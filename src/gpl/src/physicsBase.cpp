// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2018-2025, The OpenROAD Authors

// Skeleton: the coupling task fills in the private methods. The public
// surface and the checkpoint schedule are final.

#include "physicsBase.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <utility>

#include "gpl/Replace.h"
#include "nesterovBase.h"
#include "utl/Logger.h"

namespace gpl {

using utl::GPL;

PhysicsBaseVars::PhysicsBaseVars(const PlaceOptions& options)
    : overflows(options.physicsDrivenOverflows),
      thermal_weight(options.physicsThermalWeight),
      ir_weight(options.physicsIrWeight),
      reliability_weight(options.physicsReliabilityWeight),
      timing_weight(options.physicsTimingWeight),
      force_ratio(options.physicsForceRatio),
      max_inflation_ratio(options.physicsMaxInflationRatio),
      max_inflation_total(options.physicsMaxInflationTotal)
{
  model.ambient_c = options.physicsAmbientTemperature;
  model.t_ref_c = options.physicsAmbientTemperature;
  if (options.physicsSupplyVoltage > 0) {
    model.vdd_v = options.physicsSupplyVoltage;
  }
  model.pad_pitch_um = options.physicsPadPitch;
  std::ranges::sort(overflows, std::greater<int>());
}

PhysicsBase::PhysicsBase(PhysicsBaseVars vars,
                         odb::dbDatabase* db,
                         sta::dbSta* sta,
                         std::shared_ptr<NesterovBaseCommon> nbc,
                         std::vector<std::shared_ptr<NesterovBase>> nbVec,
                         utl::Logger* log)
    : vars_(std::move(vars)),
      db_(db),
      sta_(sta),
      nbc_(std::move(nbc)),
      nbVec_(std::move(nbVec)),
      log_(log)
{
  initOverflowChk();
}

PhysicsBase::~PhysicsBase() = default;

void PhysicsBase::initOverflowChk()
{
  overflow_chk_.assign(vars_.overflows.size(), false);
}

bool PhysicsBase::isPhysicsOverflow(float overflow)
{
  const int int_overflow = std::round(overflow * 100);
  if (vars_.overflows.empty() || int_overflow > vars_.overflows[0]) {
    return false;
  }
  bool need_run = false;
  for (size_t i = 0; i < vars_.overflows.size(); i++) {
    if (vars_.overflows[i] > int_overflow) {
      if (!overflow_chk_[i]) {
        overflow_chk_[i] = true;
        need_run = true;
      }
      continue;
    }
    break;
  }
  return need_run;
}

bool PhysicsBase::executePhysicsDriven(int iter,
                                       float overflow,
                                       bool apply_actuators)
{
  // Stub: coupling task implements power collection, binning, solve and
  // actuators. Keep the pass bookkeeping so the schedule is observable.
  (void) apply_actuators;
  pass_count_++;
  PhysicsSnapshot::PassRecord record;
  record.iter = iter;
  record.overflow = overflow;
  if (fields_) {
    record.summary = fields_->summary();
  }
  history_.push_back(record);
  return false;
}

void PhysicsBase::clearActuators()
{
}

void PhysicsBase::resolveSupplyVoltage()
{
}

void PhysicsBase::collectInstancePower()
{
}

void PhysicsBase::buildGrid()
{
}

void PhysicsBase::binPower(std::vector<float>& power_w) const
{
  power_w.assign(grid_.size(), 0.0f);
}

void PhysicsBase::applyFieldForce()
{
}

void PhysicsBase::applyInflation()
{
}

void PhysicsBase::applyNetWeights()
{
}

PhysicsSnapshot PhysicsBase::makeSnapshot() const
{
  PhysicsSnapshot snap;
  snap.grid = grid_;
  snap.grid_lx_dbu = grid_lx_dbu_;
  snap.grid_ly_dbu = grid_ly_dbu_;
  snap.bin_w_dbu = bin_w_dbu_;
  snap.bin_h_dbu = bin_h_dbu_;
  if (fields_) {
    snap.power_density = fields_->powerDensity();
    snap.temperature = fields_->temperature();
    snap.ir_drop = fields_->irDrop();
    snap.delay_derate = fields_->delayDerate();
    snap.stress = fields_->stress();
    snap.pad_bins = fields_->padBins();
    snap.summary = fields_->summary();
  }
  snap.passes = pass_count_;
  snap.inflated_area_percent = getInflatedAreaPercent();
  snap.history = history_;
  snap.used_sta_power = used_sta_power_;
  return snap;
}

float PhysicsBase::getInflatedAreaPercent() const
{
  if (original_movable_area_ <= 0) {
    return 0.0f;
  }
  return 100.0f * static_cast<float>(inflated_area_)
         / static_cast<float>(original_movable_area_);
}

std::string PhysicsBase::summaryString(const PhysicsSummary& s)
{
  return fmt::format(
      "power {:.4g} W, Tmax {:.1f} C (avg {:.1f}), IR drop max {:.1f} mV "
      "({:.2f}%), derate max {:.3f}, stress max {:.2f}",
      s.total_power_w,
      s.max_temperature_c,
      s.avg_temperature_c,
      1000.0f * s.max_ir_drop_v,
      s.max_ir_drop_percent,
      s.max_delay_derate,
      s.max_stress);
}

}  // namespace gpl
