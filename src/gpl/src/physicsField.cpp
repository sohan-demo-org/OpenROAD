// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2018-2025, The OpenROAD Authors

// Pure numeric coupled power / thermal / IR-drop / timing-derate /
// reliability model on a bin grid. Deterministic, single-threaded,
// no odb/sta dependency.

#include "physicsField.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>

namespace gpl {

namespace {

// Lower/upper clamps on the SOR relaxation factor. The optimal factor for a
// Poisson-like problem on an N x N grid is 2 / (1 + sin(pi / N)), which
// approaches 2 for large grids; clamping keeps the sweep stable on tiny and
// strongly screened grids.
constexpr double kMinOmega = 1.5;
constexpr double kMaxOmega = 1.98;

double optimalOmega(int cnt_x, int cnt_y)
{
  const int n = std::max(cnt_x, cnt_y);
  if (n <= 2) {
    return 1.0;
  }
  const double omega = 2.0 / (1.0 + std::sin(std::numbers::pi / n));
  return std::clamp(omega, kMinOmega, kMaxOmega);
}

// Bin indices of lattice points at pitch/2 + i * pitch along an axis of
// `cnt` bins of size `bin`; the axis center when the axis is shorter than
// a pitch.
std::vector<int> latticeBins(double pitch, double bin, int cnt)
{
  std::vector<int> bins;
  const double extent = cnt * bin;
  for (double p = pitch / 2; p < extent; p += pitch) {
    bins.push_back(std::min(static_cast<int>(p / bin), cnt - 1));
  }
  if (bins.empty() && cnt > 0) {
    bins.push_back(cnt / 2);
  }
  return bins;
}

}  // namespace

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
  user_pads_ = true;
}

bool PhysicsFields::solve()
{
  if (grid_.size() == 0 || power_w_.size() != grid_.size() || grid_.cnt_x <= 0
      || grid_.cnt_y <= 0) {
    return false;
  }
  const size_t n = grid_.size();
  power_density_.assign(n, 0.0f);
  const float area = std::max(grid_.binAreaUm2(), 1e-12f);
  for (size_t i = 0; i < n; i++) {
    power_density_[i] = power_w_[i] / area;
  }
  if (!user_pads_ || is_pad_.size() != n) {
    buildPadBins();
  }
  solveThermal();
  solveIrDrop();
  computeDerate();
  computeStress();
  computeSummary();
  cacheRanges();
  solved_ = true;
  return true;
}

void PhysicsFields::buildPadBins()
{
  const size_t n = grid_.size();
  is_pad_.assign(n, false);
  const int cx = grid_.cnt_x;
  const int cy = grid_.cnt_y;

  // Flip-chip bump lattice: pads every pad_pitch_um, starting half a pitch in
  // from the lower-left corner of the grid. An axis shorter than one pitch
  // (thin strip) gets a single row of bumps through its center.
  if (params_.pad_pitch_um > 0 && grid_.bin_w_um > 0 && grid_.bin_h_um > 0) {
    const std::vector<int> xs
        = latticeBins(params_.pad_pitch_um, grid_.bin_w_um, cx);
    const std::vector<int> ys
        = latticeBins(params_.pad_pitch_um, grid_.bin_h_um, cy);
    for (const int x : xs) {
      for (const int y : ys) {
        is_pad_[grid_.index(x, y)] = true;
      }
    }
    if (!xs.empty() && !ys.empty()) {
      return;
    }
  }

  // Wire-bond style: the outer ring of bins is the pad ring.
  for (int x = 0; x < cx; x++) {
    for (int y = 0; y < cy; y++) {
      if (x == 0 || y == 0 || x == cx - 1 || y == cy - 1) {
        is_pad_[grid_.index(x, y)] = true;
      }
    }
  }
}

// Red-black SOR sweep of the 5-point stencil
//   diag_extra * u - kx * (u[x-1] + u[x+1] - nx * u) - ky * (...) = rhs
// with Neumann (zero-flux: missing neighbors drop out) lateral edges and
// optional Dirichlet (u = 0) bins flagged in `fixed`. Returns the number of
// iterations run; `u` must be pre-sized to grid().size().
int PhysicsFields::sorSolve(const std::vector<double>& rhs,
                            double kx,
                            double ky,
                            double diag_extra,
                            const std::vector<bool>* fixed,
                            std::vector<double>& u) const
{
  const int cx = grid_.cnt_x;
  const int cy = grid_.cnt_y;
  const double omega = optimalOmega(cx, cy);
  const double tol = params_.solver_tolerance;
  const int max_iters = std::max(params_.max_solver_iters, 1);

  int iter = 0;
  for (; iter < max_iters; iter++) {
    double max_delta = 0.0;
    double max_abs = 0.0;
    for (int color = 0; color < 2; color++) {
      for (int x = 0; x < cx; x++) {
        const size_t col = grid_.index(x, 0);
        const size_t col_l = x > 0 ? grid_.index(x - 1, 0) : col;
        const size_t col_r = x + 1 < cx ? grid_.index(x + 1, 0) : col;
        const int nx = (x > 0 ? 1 : 0) + (x + 1 < cx ? 1 : 0);
        for (int y = (x + color) & 1; y < cy; y += 2) {
          const size_t idx = col + y;
          if (fixed != nullptr && (*fixed)[idx]) {
            continue;
          }
          double acc = rhs[idx];
          double diag = diag_extra + nx * kx;
          if (x > 0) {
            acc += kx * u[col_l + y];
          }
          if (x + 1 < cx) {
            acc += kx * u[col_r + y];
          }
          if (y > 0) {
            acc += ky * u[idx - 1];
            diag += ky;
          }
          if (y + 1 < cy) {
            acc += ky * u[idx + 1];
            diag += ky;
          }
          const double old = u[idx];
          const double gs = acc / diag;
          const double next = old + omega * (gs - old);
          u[idx] = next;
          max_delta = std::max(max_delta, std::abs(next - old));
          max_abs = std::max(max_abs, std::abs(next));
        }
      }
    }
    if (max_delta <= tol * std::max(max_abs, 1e-300)) {
      iter++;
      break;
    }
  }
  return iter;
}

// Steady-state screened Poisson for the temperature rise T' = T - ambient:
//   k * lap(T') - h * T' = -q
void PhysicsFields::solveThermal()
{
  const size_t n = grid_.size();
  rhs_.assign(n, 0.0);
  work_.assign(n, 0.0);
  temperature_.assign(n, params_.ambient_c);

  bool any_power = false;
  for (size_t i = 0; i < n; i++) {
    rhs_[i] = power_density_[i];
    any_power = any_power || power_density_[i] != 0.0f;
  }
  if (!any_power) {
    thermal_iters_ = 0;
    return;
  }

  const double k = std::max<double>(params_.thermal_conductivity_w_per_um_k,
                                    std::numeric_limits<double>::min());
  const double h = std::max<double>(params_.package_conductance_w_per_um2_k,
                                    std::numeric_limits<double>::min());
  const double dx = std::max<double>(grid_.bin_w_um, 1e-6);
  const double dy = std::max<double>(grid_.bin_h_um, 1e-6);
  thermal_iters_
      = sorSolve(rhs_, k / (dx * dx), k / (dy * dy), h, nullptr, work_);

  for (size_t i = 0; i < n; i++) {
    temperature_[i] = params_.ambient_c + static_cast<float>(work_[i]);
  }
}

// Resistive-sheet supply model solved for the drop D = vdd - V:
//   sigma * lap(D) = -J,  J = power density / vdd,  D = 0 at pad bins.
void PhysicsFields::solveIrDrop()
{
  const size_t n = grid_.size();
  ir_drop_.assign(n, 0.0f);
  ir_iters_ = 0;
  if (params_.vdd_v <= 0) {
    return;
  }
  const bool any_pad
      = std::find(is_pad_.begin(), is_pad_.end(), true) != is_pad_.end();
  if (!any_pad) {
    return;
  }

  rhs_.assign(n, 0.0);
  work_.assign(n, 0.0);
  bool any_power = false;
  const double vdd = params_.vdd_v;
  for (size_t i = 0; i < n; i++) {
    rhs_[i] = power_density_[i] / vdd;
    any_power = any_power || power_density_[i] != 0.0f;
  }
  if (!any_power) {
    return;
  }

  const double sigma = std::max<double>(params_.sheet_conductance_s,
                                        std::numeric_limits<double>::min());
  const double dx = std::max<double>(grid_.bin_w_um, 1e-6);
  const double dy = std::max<double>(grid_.bin_h_um, 1e-6);
  ir_iters_ = sorSolve(
      rhs_, sigma / (dx * dx), sigma / (dy * dy), 0.0, &is_pad_, work_);

  for (size_t i = 0; i < n; i++) {
    ir_drop_[i] = std::max(0.0f, static_cast<float>(work_[i]));
  }
}

void PhysicsFields::computeDerate()
{
  delay_derate_.assign(grid_.size(), 1.0f);
  for (size_t i = 0; i < grid_.size(); i++) {
    delay_derate_[i] = derateFor(temperature_[i], ir_drop_[i]);
  }
}

// S_raw = (P_density / P_ref) * arrhenius(T), normalized to max 1.
void PhysicsFields::computeStress()
{
  const size_t n = grid_.size();
  stress_.assign(n, 0.0f);
  const float p_ref
      = params_.p_ref_w_per_um2 > 0 ? params_.p_ref_w_per_um2 : 1.0f;
  float max_raw = 0.0f;
  for (size_t i = 0; i < n; i++) {
    const float raw
        = (power_density_[i] / p_ref) * arrheniusFactor(temperature_[i]);
    stress_[i] = std::max(raw, 0.0f);
    max_raw = std::max(max_raw, stress_[i]);
  }
  if (max_raw <= 0.0f) {
    std::fill(stress_.begin(), stress_.end(), 0.0f);
    return;
  }
  for (size_t i = 0; i < n; i++) {
    stress_[i] /= max_raw;
  }
}

void PhysicsFields::computeSummary()
{
  summary_ = PhysicsSummary();
  const size_t n = grid_.size();
  if (n == 0) {
    return;
  }
  double total_power = 0.0;
  double sum_t = 0.0;
  double sum_s = 0.0;
  size_t hottest = 0;
  size_t worst_ir = 0;
  size_t worst_stress = 0;
  for (size_t i = 0; i < n; i++) {
    total_power += power_w_[i];
    sum_t += temperature_[i];
    sum_s += stress_[i];
    summary_.max_power_density_w_per_um2
        = std::max(summary_.max_power_density_w_per_um2, power_density_[i]);
    summary_.max_delay_derate
        = std::max(summary_.max_delay_derate, delay_derate_[i]);
    if (temperature_[i] > temperature_[hottest]) {
      hottest = i;
    }
    if (ir_drop_[i] > ir_drop_[worst_ir]) {
      worst_ir = i;
    }
    if (stress_[i] > stress_[worst_stress]) {
      worst_stress = i;
    }
  }
  summary_.total_power_w = static_cast<float>(total_power);
  summary_.avg_temperature_c = static_cast<float>(sum_t / n);
  summary_.avg_stress = static_cast<float>(sum_s / n);
  summary_.max_temperature_c = temperature_[hottest];
  summary_.max_ir_drop_v = ir_drop_[worst_ir];
  summary_.max_stress = stress_[worst_stress];
  if (params_.vdd_v > 0) {
    summary_.max_ir_drop_percent
        = 100.0f * summary_.max_ir_drop_v / params_.vdd_v;
  }
  summary_.hottest_bin = binOf(hottest);
  summary_.worst_ir_bin = binOf(worst_ir);
  summary_.worst_stress_bin = binOf(worst_stress);
}

std::pair<int, int> PhysicsFields::binOf(size_t idx) const
{
  if (grid_.cnt_y <= 0) {
    return {0, 0};
  }
  return {static_cast<int>(idx / grid_.cnt_y),
          static_cast<int>(idx % grid_.cnt_y)};
}

PhysicsFields::FieldRange PhysicsFields::rangeOf(
    const std::vector<float>& field)
{
  FieldRange range;
  if (field.empty()) {
    return range;
  }
  const auto [min_it, max_it] = std::minmax_element(field.begin(), field.end());
  range.min = *min_it;
  range.max = *max_it;
  return range;
}

void PhysicsFields::cacheRanges()
{
  temperature_range_ = rangeOf(temperature_);
  ir_drop_range_ = rangeOf(ir_drop_);
  stress_range_ = rangeOf(stress_);
}

std::pair<float, float> PhysicsFields::gradientOf(
    const std::vector<float>& field,
    const FieldRange& range,
    int x,
    int y) const
{
  if (field.size() != grid_.size() || grid_.cnt_x == 0 || grid_.cnt_y == 0
      || x < 0 || y < 0 || x >= grid_.cnt_x || y >= grid_.cnt_y) {
    return {0.0f, 0.0f};
  }
  const float span = range.max - range.min;
  if (span <= std::numeric_limits<float>::epsilon()) {
    return {0.0f, 0.0f};
  }
  const int xl = std::max(x - 1, 0);
  const int xr = std::min(x + 1, grid_.cnt_x - 1);
  const int yl = std::max(y - 1, 0);
  const int yu = std::min(y + 1, grid_.cnt_y - 1);
  const float dx = std::max((xr - xl) * grid_.bin_w_um, 1e-6f);
  const float dy = std::max((yu - yl) * grid_.bin_h_um, 1e-6f);
  const float gx
      = (field[grid_.index(xr, y)] - field[grid_.index(xl, y)]) / (span * dx);
  const float gy
      = (field[grid_.index(x, yu)] - field[grid_.index(x, yl)]) / (span * dy);
  return {gx, gy};
}

std::pair<float, float>
PhysicsFields::gradientOf(const std::vector<float>& field, int x, int y) const
{
  return gradientOf(field, rangeOf(field), x, y);
}

std::pair<float, float> PhysicsFields::temperatureGradient(int x, int y) const
{
  return gradientOf(temperature_, temperature_range_, x, y);
}

std::pair<float, float> PhysicsFields::irDropGradient(int x, int y) const
{
  return gradientOf(ir_drop_, ir_drop_range_, x, y);
}

std::pair<float, float> PhysicsFields::stressGradient(int x, int y) const
{
  return gradientOf(stress_, stress_range_, x, y);
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
