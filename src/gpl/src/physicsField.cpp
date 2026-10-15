// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2018-2025, The OpenROAD Authors

// Skeleton: the real solvers are implemented by the field-solver task.
// This stub keeps the interface compiling and returns the trivial
// (uniform-ambient, zero-drop) solution.

#include "physicsField.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace gpl {

PhysicsFields::PhysicsFields(const PhysicsGrid& grid,
                             const PhysicsModelParams& params)
    : grid_(grid), params_(params)
{
}

void PhysicsFields::setPower(const std::vector<float>& power_w)
{
  power_w_ = power_w;
  solved_ = false;
}

void PhysicsFields::setPadBins(const std::vector<bool>& is_pad)
{
  is_pad_ = is_pad;
}

bool PhysicsFields::solve()
{
  if (grid_.size() == 0 || power_w_.size() != grid_.size()) {
    return false;
  }
  const size_t n = grid_.size();
  power_density_.assign(n, 0.0f);
  const float area = std::max(grid_.binAreaUm2(), 1e-12f);
  for (size_t i = 0; i < n; i++) {
    power_density_[i] = power_w_[i] / area;
  }
  if (is_pad_.size() != n) {
    is_pad_.assign(n, false);
  }
  solveThermal();
  solveIrDrop();
  computeDerate();
  computeStress();
  computeSummary();
  solved_ = true;
  return true;
}

void PhysicsFields::solveThermal()
{
  temperature_.assign(grid_.size(), params_.ambient_c);
}

void PhysicsFields::solveIrDrop()
{
  ir_drop_.assign(grid_.size(), 0.0f);
}

void PhysicsFields::computeDerate()
{
  delay_derate_.assign(grid_.size(), 1.0f);
  for (size_t i = 0; i < grid_.size(); i++) {
    delay_derate_[i] = derateFor(temperature_[i], ir_drop_[i]);
  }
}

void PhysicsFields::computeStress()
{
  stress_.assign(grid_.size(), 0.0f);
}

void PhysicsFields::computeSummary()
{
  summary_ = PhysicsSummary();
  for (size_t i = 0; i < grid_.size(); i++) {
    summary_.total_power_w += power_w_[i];
    summary_.max_power_density_w_per_um2
        = std::max(summary_.max_power_density_w_per_um2, power_density_[i]);
    summary_.max_temperature_c
        = std::max(summary_.max_temperature_c, temperature_[i]);
    summary_.avg_temperature_c += temperature_[i];
    summary_.max_ir_drop_v = std::max(summary_.max_ir_drop_v, ir_drop_[i]);
    summary_.max_stress = std::max(summary_.max_stress, stress_[i]);
    summary_.avg_stress += stress_[i];
    summary_.max_delay_derate
        = std::max(summary_.max_delay_derate, delay_derate_[i]);
  }
  if (grid_.size() > 0) {
    summary_.avg_temperature_c /= grid_.size();
    summary_.avg_stress /= grid_.size();
  }
  if (params_.vdd_v > 0) {
    summary_.max_ir_drop_percent
        = 100.0f * summary_.max_ir_drop_v / params_.vdd_v;
  }
}

std::pair<float, float>
PhysicsFields::gradientOf(const std::vector<float>& field, int x, int y) const
{
  if (field.size() != grid_.size() || grid_.cnt_x == 0 || grid_.cnt_y == 0) {
    return {0.0f, 0.0f};
  }
  const int xl = std::max(x - 1, 0);
  const int xr = std::min(x + 1, grid_.cnt_x - 1);
  const int yl = std::max(y - 1, 0);
  const int yu = std::min(y + 1, grid_.cnt_y - 1);
  const float dx = std::max((xr - xl) * grid_.bin_w_um, 1e-6f);
  const float dy = std::max((yu - yl) * grid_.bin_h_um, 1e-6f);
  const float gx = (normalizedAt(field, grid_.index(xr, y))
                    - normalizedAt(field, grid_.index(xl, y)))
                   / dx;
  const float gy = (normalizedAt(field, grid_.index(x, yu))
                    - normalizedAt(field, grid_.index(x, yl)))
                   / dy;
  return {gx, gy};
}

std::pair<float, float> PhysicsFields::temperatureGradient(int x, int y) const
{
  return gradientOf(temperature_, x, y);
}

std::pair<float, float> PhysicsFields::irDropGradient(int x, int y) const
{
  return gradientOf(ir_drop_, x, y);
}

std::pair<float, float> PhysicsFields::stressGradient(int x, int y) const
{
  return gradientOf(stress_, x, y);
}

float PhysicsFields::normalizedAt(const std::vector<float>& field, size_t idx)
{
  if (field.empty() || idx >= field.size()) {
    return 0.0f;
  }
  const auto [min_it, max_it] = std::minmax_element(field.begin(), field.end());
  const float range = *max_it - *min_it;
  if (range <= std::numeric_limits<float>::epsilon()) {
    return 0.0f;
  }
  return (field[idx] - *min_it) / range;
}

float PhysicsFields::derateFor(float temperature_c, float ir_drop_v) const
{
  const float thermal
      = 1.0f + params_.alpha_t_per_c * (temperature_c - params_.t_ref_c);
  const float supply
      = params_.vdd_v > 0 ? 1.0f + params_.alpha_v * (ir_drop_v / params_.vdd_v)
                          : 1.0f;
  return std::max(1.0f, thermal * supply);
}

float PhysicsFields::arrheniusFactor(float temperature_c) const
{
  constexpr float kBoltzmannEvPerK = 8.617333e-5f;
  const float t_k = temperature_c + 273.15f;
  const float t_ref_k = params_.t_ref_c + 273.15f;
  if (t_k <= 0 || t_ref_k <= 0) {
    return 1.0f;
  }
  return std::exp(-params_.activation_energy_ev / kBoltzmannEvPerK
                  * (1.0f / t_k - 1.0f / t_ref_k));
}

}  // namespace gpl
