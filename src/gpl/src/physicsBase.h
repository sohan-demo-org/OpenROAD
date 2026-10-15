// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2018-2025, The OpenROAD Authors

#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "physicsField.h"

namespace odb {
class dbDatabase;
class dbInst;
}  // namespace odb

namespace sta {
class dbSta;
}

namespace utl {
class Logger;
}

namespace gpl {

class NesterovBaseCommon;
class NesterovBase;
class GCell;
struct PlaceOptions;

// Knobs of the physics-driven mode, filled from PlaceOptions.
struct PhysicsBaseVars
{
  PhysicsBaseVars() = default;
  explicit PhysicsBaseVars(const PlaceOptions& options);

  std::vector<int> overflows{79, 64, 40, 20};
  float thermal_weight = 1.0f;
  float ir_weight = 1.0f;
  float reliability_weight = 1.0f;
  float timing_weight = 1.0f;
  float force_ratio = 0.1f;
  float max_inflation_ratio = 2.0f;
  float max_inflation_total = 0.15f;
  // Bins with stress above this are candidates for inflation.
  float inflation_stress_threshold = 0.6f;
  float inflation_coef = 2.0f;
  PhysicsModelParams model;
};

// Copy of the last solved fields, kept by Replace after the placer state is
// torn down so GUI heat maps and report_physics_placement can read it.
struct PhysicsSnapshot
{
  PhysicsGrid grid;
  // Die coordinates in DBU of the grid origin / bin size, for heat maps.
  int grid_lx_dbu = 0;
  int grid_ly_dbu = 0;
  int bin_w_dbu = 0;
  int bin_h_dbu = 0;
  std::vector<float> power_density;  // W/um^2
  std::vector<float> temperature;    // degC
  std::vector<float> ir_drop;        // V
  std::vector<float> delay_derate;   // unitless >= 1
  std::vector<float> stress;         // normalized [0, 1]
  std::vector<bool> pad_bins;
  PhysicsSummary summary;
  int passes = 0;
  float inflated_area_percent = 0;
  // Per-pass history of the headline numbers (for logs and demos).
  struct PassRecord
  {
    int iter = 0;
    float overflow = 0;
    PhysicsSummary summary;
  };
  std::vector<PassRecord> history;
};

// Couples the physics fields into Nesterov placement. One instance per
// placement run; drives three actuators at each overflow checkpoint:
//   A. per-cell field force   -> NesterovBase physics gradients
//   B. thermal/EM inflation   -> GCell::setSize on stressed cells
//   C. delay-derate net weights -> GNet::setPhysicsWeight
class PhysicsBase
{
 public:
  PhysicsBase(PhysicsBaseVars vars,
              odb::dbDatabase* db,
              sta::dbSta* sta,
              std::shared_ptr<NesterovBaseCommon> nbc,
              std::vector<std::shared_ptr<NesterovBase>> nbVec,
              utl::Logger* log);
  ~PhysicsBase();

  // Checkpoint schedule: true the first time the (unscaled) average overflow
  // drops below one of vars.overflows (percent). Mirrors
  // TimingBase::isTimingNetWeightOverflow.
  bool isPhysicsOverflow(float overflow);
  size_t getOverflowCount() const { return vars_.overflows.size(); }
  int getPassCount() const { return pass_count_; }

  // One physics pass: collect per-instance power, bin it on the top-level
  // NesterovBase bin grid, solve the fields, and (if apply_actuators) push
  // the force / inflation / net weights into the placer. Returns false when
  // nothing could be done (e.g. no movable cells).
  bool executePhysicsDriven(int iter, float overflow, bool apply_actuators);

  // Reset the actuators so a plain Nesterov run can follow (used by the final
  // reporting pass and by divergence revert).
  void clearActuators();

  // Per-instance power in watts used by the last pass (sta-based when the
  // liberty has power data, otherwise a leakage/area proxy).
  const std::vector<float>& getInstancePower() const { return inst_power_w_; }
  bool usedStaPower() const { return used_sta_power_; }

  // Supply voltage actually used (from option, SDC, liberty or default).
  float getSupplyVoltage() const { return vars_.model.vdd_v; }

  // Fields of the last pass (nullptr before the first pass).
  const PhysicsFields* getFields() const { return fields_.get(); }

  // Deep copy of the current state for Replace to keep after reset().
  PhysicsSnapshot makeSnapshot() const;

  // Total area added by inflation as a fraction of the original movable area.
  float getInflatedAreaPercent() const;

  // Human-readable summary line used by the log and the report command.
  static std::string summaryString(const PhysicsSummary& summary);

 private:
  void initOverflowChk();
  // Implemented by the coupling task.
  void resolveSupplyVoltage();
  void collectInstancePower();
  void buildGrid();
  void binPower(std::vector<float>& power_w) const;
  void applyFieldForce();
  void applyInflation();
  void applyNetWeights();

  PhysicsBaseVars vars_;
  odb::dbDatabase* db_ = nullptr;
  sta::dbSta* sta_ = nullptr;
  std::shared_ptr<NesterovBaseCommon> nbc_;
  std::vector<std::shared_ptr<NesterovBase>> nbVec_;
  utl::Logger* log_ = nullptr;

  std::vector<bool> overflow_chk_;
  int pass_count_ = 0;

  PhysicsGrid grid_;
  int grid_lx_dbu_ = 0;
  int grid_ly_dbu_ = 0;
  int bin_w_dbu_ = 0;
  int bin_h_dbu_ = 0;
  std::unique_ptr<PhysicsFields> fields_;

  // Indexed like nbc_->getGCells().
  std::vector<float> inst_power_w_;
  bool used_sta_power_ = false;

  int64_t original_movable_area_ = 0;
  int64_t inflated_area_ = 0;
  std::vector<PhysicsSnapshot::PassRecord> history_;
};

}  // namespace gpl
