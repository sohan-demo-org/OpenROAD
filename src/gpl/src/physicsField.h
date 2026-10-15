// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2018-2025, The OpenROAD Authors

#pragma once

#include <cstddef>
#include <utility>
#include <vector>

namespace gpl {

// Bin grid the physics fields live on. Flat index = x * cnt_y + y, the same
// layout the gpl FFT uses. Dimensions are in microns.
struct PhysicsGrid
{
  int cnt_x = 0;
  int cnt_y = 0;
  float bin_w_um = 0;
  float bin_h_um = 0;
  float lx_um = 0;
  float ly_um = 0;

  size_t size() const { return static_cast<size_t>(cnt_x) * cnt_y; }
  size_t index(int x, int y) const
  {
    return static_cast<size_t>(x) * cnt_y + y;
  }
  float binAreaUm2() const { return bin_w_um * bin_h_um; }
};

// Physical constants of the coupled model. Every value has a usable default
// so the mode works on any design; the Tcl options override a subset.
struct PhysicsModelParams
{
  // Thermal: k * lap(T) - h * (T - ambient) = -q
  float ambient_c = 25.0f;
  // Effective lateral thermal conductivity of the die, W/(um*K).
  // Silicon is ~1.3e-4 W/(um*K); the effective die+BEOL value is lower.
  float thermal_conductivity_w_per_um_k = 1.0e-4f;
  // Vertical heat removal to ambient through die/package, W/(um^2*K).
  float package_conductance_w_per_um2_k = 2.0e-9f;

  // Supply: sigma * lap(V) = -J, V = vdd at pads.
  float vdd_v = 1.1f;
  // Effective supply-mesh sheet conductance, S (1/ohm per square).
  float sheet_conductance_s = 2.0f;
  // Bump pitch in um; <= 0 models a boundary ring of pads (wire-bond).
  float pad_pitch_um = 0.0f;

  // Timing derate: D = (1 + alpha_t * (T - t_ref)) * (1 + alpha_v * dV/vdd)
  float t_ref_c = 25.0f;
  float alpha_t_per_c = 0.0015f;
  float alpha_v = 1.5f;

  // Reliability (Arrhenius-accelerated, current-density-driven wear-out)
  float activation_energy_ev = 0.9f;
  // Power density (W/um^2) that scores stress 1.0 at t_ref before
  // normalization.
  float p_ref_w_per_um2 = 1.0e-6f;

  // Solver controls
  int max_solver_iters = 2000;
  float solver_tolerance = 1.0e-5f;
};

// Scalar results of one solve, for logging, metrics and the report command.
struct PhysicsSummary
{
  float total_power_w = 0;
  float max_power_density_w_per_um2 = 0;
  float max_temperature_c = 0;
  float avg_temperature_c = 0;
  float max_ir_drop_v = 0;
  float max_ir_drop_percent = 0;
  float max_stress = 0;
  float avg_stress = 0;
  float max_delay_derate = 1;
  // Grid index of the hottest / worst-IR / most-stressed bins.
  std::pair<int, int> hottest_bin{0, 0};
  std::pair<int, int> worst_ir_bin{0, 0};
  std::pair<int, int> worst_stress_bin{0, 0};
};

// Pure numeric coupled power/thermal/IR/timing/reliability model on a bin
// grid. No odb/sta dependency so it is unit-testable in isolation.
//
// Usage: construct, setPower(), solve(), then read the fields. All field
// vectors have grid().size() entries indexed by grid().index(x, y).
class PhysicsFields
{
 public:
  PhysicsFields(const PhysicsGrid& grid, const PhysicsModelParams& params);

  const PhysicsGrid& grid() const { return grid_; }
  const PhysicsModelParams& params() const { return params_; }

  // Power per bin in watts (not density). Size must equal grid().size().
  void setPower(const std::vector<float>& power_w);
  // Optional: mark bins that contain a supply pad (overrides the pad model).
  void setPadBins(const std::vector<bool>& is_pad);

  // Runs thermal, IR, derate and stress in that order. Returns false if the
  // grid is empty or no power was set.
  bool solve();

  // Field accessors (valid after solve()).
  const std::vector<float>& powerDensity() const { return power_density_; }
  const std::vector<float>& temperature() const { return temperature_; }
  const std::vector<float>& irDrop() const { return ir_drop_; }
  const std::vector<float>& delayDerate() const { return delay_derate_; }
  const std::vector<float>& stress() const { return stress_; }
  const std::vector<bool>& padBins() const { return is_pad_; }

  // Central-difference gradients (per um) of the normalized fields at a bin.
  // Normalized = (field - min) / (max - min), so gradients are unit-free.
  std::pair<float, float> temperatureGradient(int x, int y) const;
  std::pair<float, float> irDropGradient(int x, int y) const;
  std::pair<float, float> stressGradient(int x, int y) const;

  const PhysicsSummary& summary() const { return summary_; }

  // Helpers exposed for tests.
  static float normalizedAt(const std::vector<float>& field, size_t idx);
  // Delay derate for a single (T, IR drop) pair under these params.
  float derateFor(float temperature_c, float ir_drop_v) const;
  // Arrhenius acceleration factor relative to t_ref.
  float arrheniusFactor(float temperature_c) const;

 private:
  void solveThermal();
  void solveIrDrop();
  void computeDerate();
  void computeStress();
  void computeSummary();
  std::pair<float, float> gradientOf(const std::vector<float>& field,
                                     int x,
                                     int y) const;

  PhysicsGrid grid_;
  PhysicsModelParams params_;

  std::vector<float> power_w_;
  std::vector<float> power_density_;
  std::vector<bool> is_pad_;
  std::vector<float> temperature_;
  std::vector<float> ir_drop_;
  std::vector<float> delay_derate_;
  std::vector<float> stress_;
  PhysicsSummary summary_;
  bool solved_ = false;
};

}  // namespace gpl
