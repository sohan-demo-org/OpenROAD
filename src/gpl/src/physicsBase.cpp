// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2018-2025, The OpenROAD Authors

// Couples the PhysicsFields solver into Nesterov placement: per-instance
// power extraction, binning on the top-level density grid, and the three
// actuators (field force, inflation, net reweighting). GPL message ids
// 200-229 are reserved for this file.

#include "physicsBase.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "db_sta/dbNetwork.hh"
#include "db_sta/dbSta.hh"
#include "gpl/Replace.h"
#include "nesterovBase.h"
#include "odb/db.h"
#include "placerBase.h"
#include "sta/Liberty.hh"
#include "sta/MinMax.hh"
#include "sta/NetworkClass.hh"
#include "sta/PowerClass.hh"
#include "sta/Scene.hh"
#include "sta/Sdc.hh"
#include "utl/Logger.h"

namespace gpl {

using utl::GPL;

namespace {

// Area proxy used when neither STA power nor liberty leakage is available.
constexpr float kAreaPowerProxyWPerUm2 = 1.0e-7f;
// Clamp of the per-cell power share that scales the field force.
constexpr float kMinPowerScale = 0.1f;
constexpr float kMaxPowerScale = 10.0f;

float finiteOrZero(float v)
{
  return std::isfinite(v) ? v : 0.0f;
}

}  // namespace

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
    vdd_from_option = true;
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
      log_(log),
      force_ratio_(vars_.force_ratio)
{
  initOverflowChk();
  pending_area_delta_.assign(nbVec_.size(), 0);
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
  if (!nbc_ || nbVec_.empty() || nbc_->getGCells().empty()) {
    return false;
  }

  if (!initialized_) {
    resolveSupplyVoltage();
    original_movable_area_ = 0;
    for (const auto& nb : nbVec_) {
      original_movable_area_ += nb->getNesterovInstsArea();
    }
    initialized_ = true;
  }

  // Power is re-collected every pass: timing-driven repair can resize or add
  // instances between checkpoints and the STA call is cheap next to a pass.
  collectInstancePower();
  buildGrid();

  std::vector<float> power_w;
  binPower(power_w);
  fields_->setPower(power_w);
  const bool solved = fields_->solve();

  if (apply_actuators) {
    pass_count_++;
  }
  PhysicsSnapshot::PassRecord record;
  record.iter = iter;
  record.overflow = overflow;
  if (solved) {
    record.summary = fields_->summary();
  }
  history_.push_back(record);

  if (!solved) {
    log_->warn(GPL,
               206,
               "Physics-driven: field solve failed on a {}x{} grid; pass "
               "skipped.",
               grid_.cnt_x,
               grid_.cnt_y);
    return false;
  }

  if (apply_actuators) {
    // Inflation first: it may remove filler GCells, and the force vector must
    // match the GCell count that updateGradients() will see.
    applyInflation();
    applyFieldForce();
    applyNetWeights();
  }

  log_->info(GPL,
             202,
             "Physics-driven {}pass {} (iter {}, overflow {:.3f}): {}",
             apply_actuators ? "" : "final ",
             apply_actuators ? pass_count_ : pass_count_ + 1,
             iter + 1,
             overflow,
             summaryString(fields_->summary()));

  if (!apply_actuators) {
    pushMetrics();
  }
  return true;
}

void PhysicsBase::pushMetrics() const
{
  const PhysicsSummary& s = fields_->summary();
  log_->metric("gpl__physics__max_temperature_c", s.max_temperature_c);
  log_->metric("gpl__physics__avg_temperature_c", s.avg_temperature_c);
  log_->metric("gpl__physics__max_ir_drop_v", s.max_ir_drop_v);
  log_->metric("gpl__physics__max_ir_drop_percent", s.max_ir_drop_percent);
  log_->metric("gpl__physics__max_stress", s.max_stress);
  log_->metric("gpl__physics__max_delay_derate", s.max_delay_derate);
  log_->metric("gpl__physics__total_power_w", s.total_power_w);
  log_->metric("gpl__physics__inflated_area_percent", getInflatedAreaPercent());
}

void PhysicsBase::clearActuators()
{
  if (!nbc_) {
    return;
  }
  for (GNet* gnet : nbc_->getGNets()) {
    gnet->setPhysicsWeight(1.0f);
  }
  nbc_->refreshDeviceNetWeights();

  // Restore the original GCell sizes and let the density model see it. The
  // reported inflated area keeps describing the last actuating pass.
  if (!tracked_sizes_.empty()) {
    const int64_t applied_inflation = inflated_area_;
    std::vector<float> no_inflation(grid_.size(), 1.0f);
    resizeTrackedGCells(no_inflation);
    syncAreasAfterResize();
    inflated_area_ = applied_inflation;
  }

  for (auto& nb : nbVec_) {
    nb->clearPhysicsGradients();
  }
}

void PhysicsBase::notifyDivergence()
{
  for (auto& nb : nbVec_) {
    nb->clearPhysicsGradients();
  }
  if (force_ratio_ <= 0) {
    return;
  }
  force_ratio_ *= 0.5f;
  log_->warn(GPL,
             205,
             "Physics-driven: divergence detected, scaling the physics "
             "force ratio down to {:.4g} for the remaining passes.",
             force_ratio_);
}

void PhysicsBase::resolveSupplyVoltage()
{
  if (vars_.vdd_from_option || sta_ == nullptr) {
    return;
  }
  sta::Scene* scene = sta_->cmdScene();
  if (scene == nullptr && !sta_->scenes().empty()) {
    scene = sta_->scenes()[0];
  }
  if (scene == nullptr) {
    return;
  }
  const sta::MinMax* max = sta::MinMax::max();
  sta::Sdc* sdc = scene->sdc();
  if (sdc != nullptr) {
    bool exists = false;
    float sdc_voltage = 0;
    sdc->voltage(max, sdc_voltage, exists);
    if (exists && sdc_voltage > 0) {
      vars_.model.vdd_v = sdc_voltage;
      return;
    }
  }
  const sta::Pvt* pvt = sdc ? sdc->operatingConditions(max) : nullptr;
  if (pvt == nullptr) {
    const sta::LibertyLibrary* lib
        = sta_->getDbNetwork()->defaultLibertyLibrary();
    if (lib != nullptr) {
      pvt = lib->defaultOperatingConditions();
    }
  }
  if (pvt != nullptr && pvt->voltage() > 0) {
    vars_.model.vdd_v = pvt->voltage();
  }
}

void PhysicsBase::collectStaPower(std::vector<float>& power, bool& any_positive)
{
  any_positive = false;
  if (sta_ == nullptr) {
    return;
  }
  sta::Scene* scene = sta_->cmdScene();
  if (scene == nullptr && !sta_->scenes().empty()) {
    scene = sta_->scenes()[0];
  }
  if (scene == nullptr) {
    return;
  }
  sta_->ensureGraph();
  sta_->ensureLevelized();
  sta::dbNetwork* network = sta_->getDbNetwork();

  const std::vector<GCell*>& gcells = nbc_->getGCells();
  for (size_t i = 0; i < gcells.size(); i++) {
    float sum = 0;
    for (const Instance* inst : gcells[i]->insts()) {
      odb::dbInst* db_inst = inst->dbInst();
      if (db_inst == nullptr) {
        continue;
      }
      sta::Instance* sta_inst = network->dbToSta(db_inst);
      if (sta_inst == nullptr || network->libertyCell(sta_inst) == nullptr) {
        continue;
      }
      const sta::PowerResult result = sta_->power(sta_inst, scene);
      sum += finiteOrZero(result.total());
    }
    power[i] = std::max(sum, 0.0f);
    any_positive = any_positive || power[i] > 0;
  }
}

void PhysicsBase::collectProxyPower(std::vector<float>& power)
{
  odb::dbBlock* block = db_->getChip()->getBlock();
  sta::dbNetwork* network = sta_ ? sta_->getDbNetwork() : nullptr;
  const std::vector<GCell*>& gcells = nbc_->getGCells();

  // First choice: liberty cell leakage.
  bool any_leakage = false;
  if (network != nullptr) {
    for (size_t i = 0; i < gcells.size(); i++) {
      float sum = 0;
      for (const Instance* inst : gcells[i]->insts()) {
        odb::dbInst* db_inst = inst->dbInst();
        if (db_inst == nullptr) {
          continue;
        }
        sta::LibertyCell* cell = network->libertyCell(db_inst);
        if (cell == nullptr) {
          continue;
        }
        float leakage = 0;
        bool exists = false;
        cell->leakagePower(leakage, exists);
        if (exists) {
          sum += std::max(finiteOrZero(leakage), 0.0f);
        }
      }
      power[i] = sum;
      any_leakage = any_leakage || sum > 0;
    }
  }
  if (any_leakage) {
    return;
  }

  // Last resort: cell area times a constant density.
  for (size_t i = 0; i < gcells.size(); i++) {
    float area_um2 = 0;
    for (const Instance* inst : gcells[i]->insts()) {
      if (inst->dbInst() == nullptr) {
        continue;
      }
      area_um2 += block->dbuAreaToMicrons(inst->getArea());
    }
    power[i] = area_um2 * kAreaPowerProxyWPerUm2;
  }
}

void PhysicsBase::collectInstancePower()
{
  const size_t count = nbc_->getGCells().size();
  inst_power_w_.assign(count, 0.0f);

  bool any_positive = false;
  collectStaPower(inst_power_w_, any_positive);
  used_sta_power_ = any_positive;
  if (!used_sta_power_) {
    std::fill(inst_power_w_.begin(), inst_power_w_.end(), 0.0f);
    collectProxyPower(inst_power_w_);
  }

  double total = 0;
  for (const float p : inst_power_w_) {
    total += p;
  }
  log_->info(GPL,
             201,
             "Physics-driven: instance power from {} for {} cells, total "
             "{:.4g} W.",
             used_sta_power_ ? "STA (liberty power tables)"
                             : "proxy (liberty leakage or cell area)",
             count,
             total);
}

void PhysicsBase::buildGrid()
{
  const NesterovBase* nb = nbVec_[0].get();
  odb::dbBlock* block = db_->getChip()->getBlock();
  const double dbu = block->getDbUnitsPerMicron();

  grid_.cnt_x = nb->getBinCntX();
  grid_.cnt_y = nb->getBinCntY();
  const std::vector<Bin>& bins = nb->getBinsConst();
  grid_lx_dbu_ = bins.empty() ? 0 : bins[0].lx();
  grid_ly_dbu_ = bins.empty() ? 0 : bins[0].ly();
  bin_w_dbu_ = static_cast<int>(std::lround(nb->getBinSizeX()));
  bin_h_dbu_ = static_cast<int>(std::lround(nb->getBinSizeY()));
  grid_.bin_w_um = static_cast<float>(nb->getBinSizeX() / dbu);
  grid_.bin_h_um = static_cast<float>(nb->getBinSizeY() / dbu);
  grid_.lx_um = static_cast<float>(grid_lx_dbu_ / dbu);
  grid_.ly_um = static_cast<float>(grid_ly_dbu_ / dbu);

  fields_ = std::make_unique<PhysicsFields>(grid_, vars_.model);
}

bool PhysicsBase::binOf(int x_dbu, int y_dbu, int& bx, int& by) const
{
  if (grid_.cnt_x <= 0 || grid_.cnt_y <= 0 || bin_w_dbu_ <= 0
      || bin_h_dbu_ <= 0) {
    return false;
  }
  bx = std::clamp((x_dbu - grid_lx_dbu_) / bin_w_dbu_, 0, grid_.cnt_x - 1);
  by = std::clamp((y_dbu - grid_ly_dbu_) / bin_h_dbu_, 0, grid_.cnt_y - 1);
  return true;
}

void PhysicsBase::binPower(std::vector<float>& power_w) const
{
  power_w.assign(grid_.size(), 0.0f);
  const std::vector<GCell*>& gcells = nbc_->getGCells();
  const size_t count = std::min(gcells.size(), inst_power_w_.size());
  for (size_t i = 0; i < count; i++) {
    const float p = inst_power_w_[i];
    if (p <= 0) {
      continue;
    }
    int bx = 0;
    int by = 0;
    if (!binOf(gcells[i]->dCx(), gcells[i]->dCy(), bx, by)) {
      continue;
    }
    power_w[grid_.index(bx, by)] += p;
  }
}

float PhysicsBase::cellPowerScale(const GCell* gcell, float mean_power) const
{
  if (mean_power <= 0) {
    return 1.0f;
  }
  const size_t idx = nbc_->getGCellIndex(gcell);
  if (idx >= inst_power_w_.size()) {
    return 1.0f;
  }
  return std::clamp(
      inst_power_w_[idx] / mean_power, kMinPowerScale, kMaxPowerScale);
}

void PhysicsBase::applyFieldForce()
{
  if (force_ratio_ <= 0) {
    for (auto& nb : nbVec_) {
      nb->clearPhysicsGradients();
    }
    return;
  }

  // Mean power over powered std cells, for the per-cell power share.
  double power_sum = 0;
  size_t powered = 0;
  const std::vector<GCell*>& gcells = nbc_->getGCells();
  for (size_t i = 0; i < std::min(gcells.size(), inst_power_w_.size()); i++) {
    if (gcells[i]->isStdInstance() && inst_power_w_[i] > 0) {
      power_sum += inst_power_w_[i];
      powered++;
    }
  }
  const float mean_power
      = powered > 0 ? static_cast<float>(power_sum / powered) : 0.0f;

  for (auto& nb : nbVec_) {
    const std::vector<GCellHandle>& handles = nb->getGCells();
    std::vector<FloatPoint> force(handles.size(), FloatPoint(0, 0));
    double force_sum = 0;
    for (size_t i = 0; i < handles.size(); i++) {
      const GCell* gcell = handles[i];
      if (!gcell->isStdInstance() || gcell->isLocked()) {
        continue;
      }
      int bx = 0;
      int by = 0;
      if (!binOf(gcell->dCx(), gcell->dCy(), bx, by)) {
        continue;
      }
      const auto [gtx, gty] = fields_->temperatureGradient(bx, by);
      const auto [gvx, gvy] = fields_->irDropGradient(bx, by);
      const auto [gsx, gsy] = fields_->stressGradient(bx, by);
      const float scale = cellPowerScale(gcell, mean_power);
      float fx = -(vars_.thermal_weight * gtx + vars_.ir_weight * gvx
                   + vars_.reliability_weight * gsx)
                 * scale;
      float fy = -(vars_.thermal_weight * gty + vars_.ir_weight * gvy
                   + vars_.reliability_weight * gsy)
                 * scale;
      fx = finiteOrZero(fx);
      fy = finiteOrZero(fy);
      force[i] = FloatPoint(fx, fy);
      force_sum += std::fabs(fx) + std::fabs(fy);
    }

    const float wl_grad_sum = nb->getWireLengthGradSum();
    if (force_sum <= 0 || !(wl_grad_sum > 0) || !std::isfinite(wl_grad_sum)) {
      nb->clearPhysicsGradients();
      continue;
    }
    // Normalize so that sum|coef * force| == force_ratio * sum|wl grad|.
    const float coef
        = static_cast<float>(force_ratio_ * wl_grad_sum / force_sum);
    nb->setPhysicsGradients(std::move(force), finiteOrZero(coef));
  }
}

void PhysicsBase::resizeTrackedGCells(const std::vector<float>& bin_ratio)
{
  inflated_area_ = 0;
  std::fill(pending_area_delta_.begin(), pending_area_delta_.end(), 0);

  for (size_t nb_index = 0; nb_index < nbVec_.size(); nb_index++) {
    for (const GCellHandle& handle : nbVec_[nb_index]->getGCells()) {
      if (!handle->isStdInstance() || !handle.isNesterovBaseCommon()) {
        continue;
      }
      GCell* gcell = nbc_->getGCellByIndex(handle.getStorageIndex());
      odb::dbInst* key = gcell->insts()[0]->dbInst();
      if (key == nullptr) {
        continue;
      }

      float ratio = 1.0f;
      int bx = 0;
      int by = 0;
      if (!bin_ratio.empty() && binOf(gcell->dCx(), gcell->dCy(), bx, by)) {
        ratio = bin_ratio[grid_.index(bx, by)];
      }

      auto it = tracked_sizes_.find(key);
      if (it == tracked_sizes_.end()) {
        if (ratio <= 1.0f) {
          continue;
        }
        it = tracked_sizes_
                 .emplace(
                     key,
                     TrackedSize{
                         gcell->dx(), gcell->dy(), gcell->dx(), gcell->dy()})
                 .first;
      } else if (gcell->dx() != it->second.set_dx
                 || gcell->dy() != it->second.set_dy) {
        // Someone else (timing-driven repair) resized the cell since our
        // last pass: treat the new size as the original.
        it->second.orig_dx = gcell->dx();
        it->second.orig_dy = gcell->dy();
      }
      TrackedSize& tracked = it->second;

      const float side = std::sqrt(std::max(ratio, 1.0f));
      const int new_dx = static_cast<int>(std::round(tracked.orig_dx * side));
      const int new_dy = static_cast<int>(std::round(tracked.orig_dy * side));
      const int64_t prev_area = static_cast<int64_t>(gcell->dx()) * gcell->dy();
      const int64_t orig_area
          = static_cast<int64_t>(tracked.orig_dx) * tracked.orig_dy;
      if (new_dx != gcell->dx() || new_dy != gcell->dy()) {
        gcell->setSize(new_dx, new_dy, GCell::GCellChange::kRoutability);
      }
      tracked.set_dx = gcell->dx();
      tracked.set_dy = gcell->dy();
      const int64_t new_area = static_cast<int64_t>(gcell->dx()) * gcell->dy();
      pending_area_delta_[nb_index] += new_area - prev_area;
      inflated_area_ += new_area - orig_area;
    }
  }
}

void PhysicsBase::syncAreasAfterResize()
{
  // Same bookkeeping as RouteBase::routability / the timing-driven path:
  // absorb added area by removing fillers, hand any remainder (or a net
  // shrink) to the target density, then refresh areas and density sizes.
  for (size_t i = 0; i < nbVec_.size(); i++) {
    NesterovBase* nb = nbVec_[i].get();
    const int64_t delta = pending_area_delta_[i];
    if (delta == 0) {
      continue;
    }
    if (delta > 0) {
      nb->cutFillerCells(delta);
    } else {
      nb->updateAreas();
      const float white_space = nb->getWhiteSpaceArea();
      if (white_space > 0) {
        nb->setTargetDensity(static_cast<float>(nb->getNesterovInstsArea()
                                                + nb->getTotalFillerArea())
                             / white_space);
        nb->setMovableArea(nb->getTargetDensity() * white_space);
      }
    }
    nb->updateAreas();
    nb->updateDensitySize();
    nb->checkConsistency();
    pending_area_delta_[i] = 0;
  }
}

void PhysicsBase::applyInflation()
{
  const bool budget_off = vars_.max_inflation_total <= 0
                          || vars_.max_inflation_ratio <= 1.0f
                          || original_movable_area_ <= 0;
  std::vector<float> bin_ratio(grid_.size(), 1.0f);
  int stressed_bins = 0;
  if (!budget_off) {
    const std::vector<float>& stress = fields_->stress();
    for (size_t i = 0; i < grid_.size() && i < stress.size(); i++) {
      const float s = stress[i];
      if (!std::isfinite(s) || s <= vars_.inflation_stress_threshold) {
        continue;
      }
      bin_ratio[i] = std::min(
          vars_.max_inflation_ratio,
          1.0f + vars_.inflation_coef * (s - vars_.inflation_stress_threshold));
      stressed_bins++;
    }
  }
  if (stressed_bins == 0 && tracked_sizes_.empty()) {
    return;
  }

  // Scale ratios down uniformly if the wanted area exceeds the budget
  // (RouteBase::scaleInflationToBudget).
  if (stressed_bins > 0) {
    double wanted = 0;
    for (const auto& nb : nbVec_) {
      for (const GCellHandle& handle : nb->getGCells()) {
        if (!handle->isStdInstance()) {
          continue;
        }
        int bx = 0;
        int by = 0;
        if (!binOf(handle->dCx(), handle->dCy(), bx, by)) {
          continue;
        }
        const float ratio = bin_ratio[grid_.index(bx, by)];
        if (ratio <= 1.0f) {
          continue;
        }
        const GCell* gcell = handle;
        auto it = tracked_sizes_.find(gcell->insts()[0]->dbInst());
        const double area
            = it != tracked_sizes_.end()
                  ? static_cast<double>(it->second.orig_dx) * it->second.orig_dy
                  : static_cast<double>(gcell->dx()) * gcell->dy();
        wanted += area * (ratio - 1.0);
      }
    }
    const double budget = static_cast<double>(original_movable_area_)
                          * vars_.max_inflation_total;
    if (wanted > budget && wanted > 0) {
      const float scale = static_cast<float>(budget / wanted);
      for (float& ratio : bin_ratio) {
        if (ratio > 1.0f) {
          ratio = 1.0f + scale * (ratio - 1.0f);
        }
      }
    }
  }

  resizeTrackedGCells(bin_ratio);
  syncAreasAfterResize();

  int inflated_count = 0;
  for (const auto& [inst, tracked] : tracked_sizes_) {
    if (tracked.set_dx > tracked.orig_dx || tracked.set_dy > tracked.orig_dy) {
      inflated_count++;
    }
  }
  log_->info(GPL,
             203,
             "Physics-driven: {} stressed bins, {} cells inflated, inflated "
             "area {:.2f}% of movable area (budget {:.2f}%).",
             stressed_bins,
             inflated_count,
             getInflatedAreaPercent(),
             100.0f * vars_.max_inflation_total);
}

void PhysicsBase::applyNetWeights()
{
  const std::vector<float>& derate = fields_->delayDerate();
  for (GNet* gnet : nbc_->getGNets()) {
    float weight = 1.0f;
    if (vars_.timing_weight > 0 && gnet->getGPins().size() > 1
        && derate.size() == grid_.size()) {
      float max_derate = 1.0f;
      for (const GPin* gpin : gnet->getGPins()) {
        const GCell* gcell = gpin->getGCell();
        if (gcell == nullptr || !gcell->isInstance()) {
          continue;
        }
        int bx = 0;
        int by = 0;
        if (!binOf(gcell->dCx(), gcell->dCy(), bx, by)) {
          continue;
        }
        max_derate = std::max(max_derate, derate[grid_.index(bx, by)]);
      }
      weight = std::pow(max_derate, vars_.timing_weight);
      if (!std::isfinite(weight) || weight < 1.0f) {
        weight = 1.0f;
      }
    }
    gnet->setPhysicsWeight(weight);
  }
  nbc_->refreshDeviceNetWeights();
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
      "power {:.4g} W (max bin density {:.3g} W/um^2), Tmax {:.1f} C "
      "(avg {:.1f}), IR drop max {:.1f} mV ({:.2f}%), derate max {:.3f}, "
      "stress max {:.2f}",
      s.total_power_w,
      s.max_power_density_w_per_um2,
      s.max_temperature_c,
      s.avg_temperature_c,
      1000.0f * s.max_ir_drop_v,
      s.max_ir_drop_percent,
      s.max_delay_derate,
      s.max_stress);
}

}  // namespace gpl
