// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2018-2025, The OpenROAD Authors

#include <chrono>
#include <cmath>
#include <cstdio>
#include <vector>

#include "gtest/gtest.h"
#include "src/gpl/src/physicsField.h"

namespace gpl {
namespace {

PhysicsGrid makeGrid(int cnt_x,
                     int cnt_y,
                     float bin_w = 10.0f,
                     float bin_h = 10.0f)
{
  PhysicsGrid grid;
  grid.cnt_x = cnt_x;
  grid.cnt_y = cnt_y;
  grid.bin_w_um = bin_w;
  grid.bin_h_um = bin_h;
  grid.lx_um = 0;
  grid.ly_um = 0;
  return grid;
}

PhysicsModelParams tightParams()
{
  PhysicsModelParams params;
  params.solver_tolerance = 1e-7f;
  params.max_solver_iters = 20000;
  return params;
}

void expectZeroPowerInvariants(const PhysicsFields& fields)
{
  const auto& params = fields.params();
  for (size_t i = 0; i < fields.grid().size(); i++) {
    EXPECT_EQ(fields.temperature()[i], params.ambient_c);
    EXPECT_EQ(fields.irDrop()[i], 0.0f);
    EXPECT_EQ(fields.delayDerate()[i], 1.0f);
    EXPECT_EQ(fields.stress()[i], 0.0f);
  }
  const PhysicsSummary& s = fields.summary();
  EXPECT_EQ(s.total_power_w, 0.0f);
  EXPECT_EQ(s.max_power_density_w_per_um2, 0.0f);
  EXPECT_EQ(s.max_temperature_c, params.ambient_c);
  EXPECT_EQ(s.avg_temperature_c, params.ambient_c);
  EXPECT_EQ(s.max_ir_drop_v, 0.0f);
  EXPECT_EQ(s.max_ir_drop_percent, 0.0f);
  EXPECT_EQ(s.max_stress, 0.0f);
  EXPECT_EQ(s.avg_stress, 0.0f);
  EXPECT_EQ(s.max_delay_derate, 1.0f);
}

TEST(PhysicsFieldTest, EmptyOrMismatchedInputFails)
{
  PhysicsFields empty(makeGrid(0, 0), PhysicsModelParams());
  EXPECT_FALSE(empty.solve());

  PhysicsFields fields(makeGrid(4, 4), PhysicsModelParams());
  EXPECT_FALSE(fields.solve());
  fields.setPower(std::vector<float>(3, 0.0f));
  EXPECT_FALSE(fields.solve());
  fields.setPower(std::vector<float>(16, 0.0f));
  EXPECT_TRUE(fields.solve());
  EXPECT_TRUE(fields.solved());
}

TEST(PhysicsFieldTest, ZeroPowerGivesTrivialFields)
{
  const PhysicsGrid grid = makeGrid(16, 12);
  PhysicsFields fields(grid, PhysicsModelParams());
  fields.setPower(std::vector<float>(grid.size(), 0.0f));
  ASSERT_TRUE(fields.solve());
  expectZeroPowerInvariants(fields);
  EXPECT_EQ(fields.thermalIterations(), 0);
  EXPECT_EQ(fields.irDropIterations(), 0);
  EXPECT_EQ(fields.temperatureGradient(5, 5), std::make_pair(0.0f, 0.0f));
  EXPECT_EQ(fields.irDropGradient(5, 5), std::make_pair(0.0f, 0.0f));
  EXPECT_EQ(fields.stressGradient(5, 5), std::make_pair(0.0f, 0.0f));
}

TEST(PhysicsFieldTest, SingleHotBinThermalProfile)
{
  const int n = 33;
  const int c = n / 2;
  const PhysicsGrid grid = makeGrid(n, n);
  PhysicsFields fields(grid, tightParams());
  std::vector<float> power(grid.size(), 0.0f);
  const float hot_w = 2.5e-3f;
  power[grid.index(c, c)] = hot_w;
  fields.setPower(power);
  ASSERT_TRUE(fields.solve());
  EXPECT_GT(fields.thermalIterations(), 0);

  const auto& t = fields.temperature();
  const float t_max = t[grid.index(c, c)];
  EXPECT_GT(t_max, fields.params().ambient_c);
  for (size_t i = 0; i < grid.size(); i++) {
    EXPECT_LE(t[i], t_max);
    EXPECT_GE(t[i], fields.params().ambient_c);
  }
  // Monotonically non-increasing with distance along the center row.
  const float slack = 1e-6f * (t_max - fields.params().ambient_c);
  for (int x = c; x + 1 < n; x++) {
    EXPECT_GE(t[grid.index(x, c)] + slack, t[grid.index(x + 1, c)])
        << "x=" << x;
  }
  // Left/right and up/down symmetry.
  for (int d = 1; d <= c; d++) {
    EXPECT_NEAR(t[grid.index(c - d, c)], t[grid.index(c + d, c)], slack)
        << "d=" << d;
    EXPECT_NEAR(t[grid.index(c, c - d)], t[grid.index(c, c + d)], slack)
        << "d=" << d;
  }

  const PhysicsSummary& s = fields.summary();
  EXPECT_NEAR(s.total_power_w, hot_w, 1e-9f);
  EXPECT_EQ(s.hottest_bin, std::make_pair(c, c));
  EXPECT_EQ(s.worst_stress_bin, std::make_pair(c, c));
  EXPECT_FLOAT_EQ(s.max_temperature_c, t_max);
  EXPECT_FLOAT_EQ(s.max_power_density_w_per_um2, hot_w / grid.binAreaUm2());
  EXPECT_GT(s.max_delay_derate, 1.0f);
  EXPECT_FLOAT_EQ(s.max_stress, 1.0f);
}

TEST(PhysicsFieldTest, IrDropBoundaryPadsVersusBumpLattice)
{
  const int n = 21;
  const PhysicsGrid grid = makeGrid(n, n);
  std::vector<float> power(grid.size(), 1e-5f);

  PhysicsModelParams ring = tightParams();
  ring.pad_pitch_um = 0;
  PhysicsFields ring_fields(grid, ring);
  ring_fields.setPower(power);
  ASSERT_TRUE(ring_fields.solve());
  EXPECT_GT(ring_fields.irDropIterations(), 0);

  const auto& drop = ring_fields.irDrop();
  const auto& pads = ring_fields.padBins();
  float interior_max = 0;
  for (int x = 0; x < n; x++) {
    for (int y = 0; y < n; y++) {
      const size_t idx = grid.index(x, y);
      const bool on_ring = x == 0 || y == 0 || x == n - 1 || y == n - 1;
      EXPECT_EQ(pads[idx], on_ring);
      if (on_ring) {
        EXPECT_EQ(drop[idx], 0.0f);
      } else {
        EXPECT_GT(drop[idx], 0.0f);
        interior_max = std::max(interior_max, drop[idx]);
      }
    }
  }
  const PhysicsSummary& s = ring_fields.summary();
  EXPECT_FLOAT_EQ(s.max_ir_drop_v, interior_max);
  EXPECT_EQ(s.worst_ir_bin, std::make_pair(n / 2, n / 2));
  EXPECT_FLOAT_EQ(s.max_ir_drop_percent, 100.0f * interior_max / ring.vdd_v);

  PhysicsModelParams bumps = tightParams();
  bumps.pad_pitch_um = 2 * grid.bin_w_um;
  PhysicsFields bump_fields(grid, bumps);
  bump_fields.setPower(power);
  ASSERT_TRUE(bump_fields.solve());
  // Pads at the odd bins in both directions.
  EXPECT_TRUE(bump_fields.padBins()[grid.index(1, 1)]);
  EXPECT_TRUE(bump_fields.padBins()[grid.index(3, 1)]);
  EXPECT_FALSE(bump_fields.padBins()[grid.index(0, 0)]);
  EXPECT_FALSE(bump_fields.padBins()[grid.index(2, 2)]);
  EXPECT_LT(bump_fields.summary().max_ir_drop_v, 0.1f * interior_max);
  EXPECT_GT(bump_fields.summary().max_ir_drop_v, 0.0f);
}

TEST(PhysicsFieldTest, UserPadBinsOverridePadModel)
{
  const PhysicsGrid grid = makeGrid(9, 9);
  PhysicsFields fields(grid, tightParams());
  std::vector<bool> pads(grid.size(), false);
  pads[grid.index(0, 0)] = true;
  fields.setPadBins(pads);
  fields.setPower(std::vector<float>(grid.size(), 1e-5f));
  ASSERT_TRUE(fields.solve());
  EXPECT_EQ(fields.padBins(), pads);
  EXPECT_EQ(fields.irDrop()[grid.index(0, 0)], 0.0f);
  EXPECT_EQ(fields.summary().worst_ir_bin, std::make_pair(8, 8));
}

TEST(PhysicsFieldTest, NoPadBinsGivesZeroDrop)
{
  const PhysicsGrid grid = makeGrid(5, 5);
  PhysicsFields fields(grid, tightParams());
  fields.setPadBins(std::vector<bool>(grid.size(), false));
  fields.setPower(std::vector<float>(grid.size(), 1e-5f));
  ASSERT_TRUE(fields.solve());
  for (float d : fields.irDrop()) {
    EXPECT_EQ(d, 0.0f);
  }
  EXPECT_GT(fields.summary().max_temperature_c, fields.params().ambient_c);
}

TEST(PhysicsFieldTest, DoublingPowerDoublesRiseAndDrop)
{
  const PhysicsGrid grid = makeGrid(24, 18);
  std::vector<float> power(grid.size());
  for (size_t i = 0; i < grid.size(); i++) {
    power[i] = 1e-6f * static_cast<float>((i * 7919) % 13);
  }
  std::vector<float> doubled(power);
  for (float& p : doubled) {
    p *= 2.0f;
  }

  PhysicsFields a(grid, tightParams());
  a.setPower(power);
  ASSERT_TRUE(a.solve());
  PhysicsFields b(grid, tightParams());
  b.setPower(doubled);
  ASSERT_TRUE(b.solve());

  const float ambient = a.params().ambient_c;
  for (size_t i = 0; i < grid.size(); i++) {
    const float rise_a = a.temperature()[i] - ambient;
    const float rise_b = b.temperature()[i] - ambient;
    EXPECT_NEAR(rise_b, 2.0f * rise_a, 1e-3f * std::abs(rise_b) + 1e-6f);
    EXPECT_NEAR(
        b.irDrop()[i], 2.0f * a.irDrop()[i], 1e-3f * b.irDrop()[i] + 1e-9f);
  }
  EXPECT_NEAR(b.summary().total_power_w,
              2.0f * a.summary().total_power_w,
              1e-6f * b.summary().total_power_w);
}

TEST(PhysicsFieldTest, DerateAndArrheniusSpotValues)
{
  PhysicsModelParams params;
  params.t_ref_c = 25.0f;
  params.alpha_t_per_c = 0.0015f;
  params.alpha_v = 1.5f;
  params.vdd_v = 1.0f;
  params.activation_energy_ev = 0.9f;
  PhysicsFields fields(makeGrid(1, 1), params);

  EXPECT_FLOAT_EQ(fields.derateFor(25.0f, 0.0f), 1.0f);
  EXPECT_FLOAT_EQ(fields.derateFor(125.0f, 0.0f), 1.15f);
  EXPECT_FLOAT_EQ(fields.derateFor(25.0f, 0.1f), 1.15f);
  EXPECT_NEAR(fields.derateFor(125.0f, 0.1f), 1.15f * 1.15f, 1e-6f);
  // Below reference temperature the derate is clamped to 1.
  EXPECT_FLOAT_EQ(fields.derateFor(-40.0f, 0.0f), 1.0f);

  EXPECT_FLOAT_EQ(fields.arrheniusFactor(25.0f), 1.0f);
  // exp(-Ea/kB * (1/398.15 - 1/298.15)) for Ea = 0.9 eV.
  const double expected
      = std::exp(-0.9 / 8.617333e-5 * (1.0 / 398.15 - 1.0 / 298.15));
  EXPECT_NEAR(fields.arrheniusFactor(125.0f), expected, 1e-3 * expected);
  EXPECT_GT(fields.arrheniusFactor(125.0f), 1000.0f);
  EXPECT_LT(fields.arrheniusFactor(0.0f), 1.0f);
}

TEST(PhysicsFieldTest, StressNormalizedToMaxOne)
{
  const PhysicsGrid grid = makeGrid(8, 8);
  PhysicsFields fields(grid, tightParams());
  std::vector<float> power(grid.size(), 1e-6f);
  power[grid.index(2, 5)] = 4e-5f;
  fields.setPower(power);
  ASSERT_TRUE(fields.solve());
  float max_stress = 0;
  for (float s : fields.stress()) {
    EXPECT_GE(s, 0.0f);
    EXPECT_LE(s, 1.0f);
    max_stress = std::max(max_stress, s);
  }
  EXPECT_FLOAT_EQ(max_stress, 1.0f);
  EXPECT_FLOAT_EQ(fields.stress()[grid.index(2, 5)], 1.0f);
  EXPECT_EQ(fields.summary().worst_stress_bin, std::make_pair(2, 5));
  EXPECT_FLOAT_EQ(fields.summary().max_stress, 1.0f);
}

TEST(PhysicsFieldTest, GradientPointsAwayFromHotSpot)
{
  const int n = 15;
  const int c = n / 2;
  const PhysicsGrid grid = makeGrid(n, n);
  PhysicsFields fields(grid, tightParams());
  std::vector<float> power(grid.size(), 0.0f);
  power[grid.index(c, c)] = 1e-3f;
  fields.setPower(power);
  ASSERT_TRUE(fields.solve());

  EXPECT_LT(fields.temperatureGradient(c + 2, c).first, 0.0f);
  EXPECT_GT(fields.temperatureGradient(c - 2, c).first, 0.0f);
  EXPECT_LT(fields.temperatureGradient(c, c + 2).second, 0.0f);
  EXPECT_NEAR(fields.temperatureGradient(c, c).first, 0.0f, 1e-6f);
  EXPECT_LT(fields.stressGradient(c + 1, c).first, 0.0f);
  // Hot spot is the worst IR bin too (ring pads): drop falls to the right.
  EXPECT_LT(fields.irDropGradient(c + 2, c).first, 0.0f);

  // Cached ranges match the static helper.
  const auto& t = fields.temperature();
  const float expected_gx
      = (PhysicsFields::normalizedAt(t, grid.index(c + 3, c))
         - PhysicsFields::normalizedAt(t, grid.index(c + 1, c)))
        / (2 * grid.bin_w_um);
  EXPECT_NEAR(fields.temperatureGradient(c + 2, c).first,
              expected_gx,
              1e-6f * std::abs(expected_gx));
  // Out-of-range queries are safe.
  EXPECT_EQ(fields.temperatureGradient(-1, 0), std::make_pair(0.0f, 0.0f));
  EXPECT_EQ(fields.temperatureGradient(n, 0), std::make_pair(0.0f, 0.0f));
}

TEST(PhysicsFieldTest, NonSquareBinsAndThinGrids)
{
  {
    const PhysicsGrid grid = makeGrid(12, 7, 3.0f, 11.0f);
    PhysicsFields fields(grid, tightParams());
    fields.setPower(std::vector<float>(grid.size(), 0.0f));
    ASSERT_TRUE(fields.solve());
    expectZeroPowerInvariants(fields);

    std::vector<float> power(grid.size(), 0.0f);
    power[grid.index(6, 3)] = 1e-4f;
    fields.setPower(power);
    ASSERT_TRUE(fields.solve());
    EXPECT_EQ(fields.summary().hottest_bin, std::make_pair(6, 3));
    EXPECT_GT(fields.summary().max_ir_drop_v, 0.0f);
    EXPECT_LT(fields.temperatureGradient(8, 3).first, 0.0f);
  }
  for (const PhysicsGrid& grid : {makeGrid(1, 40, 5.0f, 2.0f),
                                  makeGrid(40, 1, 2.0f, 5.0f),
                                  makeGrid(1, 1)}) {
    PhysicsFields fields(grid, tightParams());
    fields.setPower(std::vector<float>(grid.size(), 0.0f));
    ASSERT_TRUE(fields.solve());
    expectZeroPowerInvariants(fields);

    std::vector<float> power(grid.size(), 2e-6f);
    fields.setPower(power);
    ASSERT_TRUE(fields.solve());
    for (size_t i = 0; i < grid.size(); i++) {
      EXPECT_TRUE(std::isfinite(fields.temperature()[i]));
      EXPECT_TRUE(std::isfinite(fields.irDrop()[i]));
      EXPECT_GE(fields.irDrop()[i], 0.0f);
      EXPECT_GT(fields.temperature()[i], fields.params().ambient_c);
    }
    EXPECT_NEAR(fields.summary().total_power_w, 2e-6f * grid.size(), 1e-9f);
  }
  // A 1xN strip with ring pads has every bin on the ring: no drop anywhere.
  {
    const PhysicsGrid grid = makeGrid(1, 10);
    PhysicsFields fields(grid, tightParams());
    fields.setPower(std::vector<float>(grid.size(), 1e-6f));
    ASSERT_TRUE(fields.solve());
    EXPECT_EQ(fields.summary().max_ir_drop_v, 0.0f);
  }
  // Bump lattice on a 1xN strip with a pitch of 4 bins.
  {
    const PhysicsGrid grid = makeGrid(1, 16, 10.0f, 10.0f);
    PhysicsModelParams params = tightParams();
    params.pad_pitch_um = 40.0f;
    PhysicsFields fields(grid, params);
    fields.setPower(std::vector<float>(grid.size(), 1e-6f));
    ASSERT_TRUE(fields.solve());
    EXPECT_TRUE(fields.padBins()[grid.index(0, 2)]);
    EXPECT_TRUE(fields.padBins()[grid.index(0, 6)]);
    EXPECT_FALSE(fields.padBins()[grid.index(0, 0)]);
    EXPECT_GT(fields.summary().max_ir_drop_v, 0.0f);
    EXPECT_EQ(fields.irDrop()[grid.index(0, 2)], 0.0f);
  }
}

// Discrete residual of the 5-point stencil with Neumann edges:
//   diag * u - kx * (u[x-1] + u[x+1]) - ky * (u[y-1] + u[y+1]) - rhs
// returned as max |residual| / max |rhs|.
double stencilResidual(const PhysicsGrid& grid,
                       const std::vector<float>& u,
                       const std::vector<float>& rhs,
                       double kx,
                       double ky,
                       double diag_extra,
                       const std::vector<bool>* fixed)
{
  double max_res = 0;
  double max_rhs = 0;
  for (int x = 0; x < grid.cnt_x; x++) {
    for (int y = 0; y < grid.cnt_y; y++) {
      const size_t idx = grid.index(x, y);
      max_rhs = std::max<double>(max_rhs, std::abs(rhs[idx]));
      if (fixed != nullptr && (*fixed)[idx]) {
        continue;
      }
      double acc = diag_extra * u[idx];
      if (x > 0) {
        acc += kx * (u[idx] - u[grid.index(x - 1, y)]);
      }
      if (x + 1 < grid.cnt_x) {
        acc += kx * (u[idx] - u[grid.index(x + 1, y)]);
      }
      if (y > 0) {
        acc += ky * (u[idx] - u[grid.index(x, y - 1)]);
      }
      if (y + 1 < grid.cnt_y) {
        acc += ky * (u[idx] - u[grid.index(x, y + 1)]);
      }
      max_res = std::max(max_res, std::abs(acc - rhs[idx]));
    }
  }
  return max_res / std::max(max_rhs, 1e-300);
}

TEST(PhysicsFieldTest, SolutionsSatisfyDiscretePdes)
{
  const PhysicsGrid grid = makeGrid(40, 28, 6.0f, 4.0f);
  PhysicsModelParams params = tightParams();
  params.pad_pitch_um = 48.0f;
  PhysicsFields fields(grid, params);
  std::vector<float> power(grid.size());
  for (size_t i = 0; i < grid.size(); i++) {
    power[i] = 1e-6f * static_cast<float>(1 + (i * 31) % 17);
  }
  fields.setPower(power);
  ASSERT_TRUE(fields.solve());

  const double dx2 = grid.bin_w_um * grid.bin_w_um;
  const double dy2 = grid.bin_h_um * grid.bin_h_um;

  std::vector<float> rise(grid.size());
  for (size_t i = 0; i < grid.size(); i++) {
    rise[i] = fields.temperature()[i] - params.ambient_c;
  }
  const double k = params.thermal_conductivity_w_per_um_k;
  EXPECT_LT(stencilResidual(grid,
                            rise,
                            fields.powerDensity(),
                            k / dx2,
                            k / dy2,
                            params.package_conductance_w_per_um2_k,
                            nullptr),
            1e-3);

  std::vector<float> current(grid.size());
  for (size_t i = 0; i < grid.size(); i++) {
    current[i] = fields.powerDensity()[i] / params.vdd_v;
  }
  const double sigma = params.sheet_conductance_s;
  EXPECT_LT(stencilResidual(grid,
                            fields.irDrop(),
                            current,
                            sigma / dx2,
                            sigma / dy2,
                            0.0,
                            &fields.padBins()),
            1e-3);
}

TEST(PhysicsFieldTest, LargeGridSolvesQuickly)
{
  const int n = 256;
  const PhysicsGrid grid = makeGrid(n, n, 4.0f, 4.0f);
  PhysicsFields fields(grid, PhysicsModelParams());
  std::vector<float> power(grid.size());
  for (int x = 0; x < n; x++) {
    for (int y = 0; y < n; y++) {
      const float dx = (x - n / 2) / 40.0f;
      const float dy = (y - n / 3) / 40.0f;
      power[grid.index(x, y)]
          = 1e-7f * (1.0f + 20.0f * std::exp(-(dx * dx + dy * dy)));
    }
  }
  fields.setPower(power);
  const auto start = std::chrono::steady_clock::now();
  ASSERT_TRUE(fields.solve());
  const double secs
      = std::chrono::duration<double>(std::chrono::steady_clock::now() - start)
            .count();
  std::printf("256x256 solve: %.3f s (thermal %d iters, IR %d iters)\n",
              secs,
              fields.thermalIterations(),
              fields.irDropIterations());
  EXPECT_LT(fields.thermalIterations(), fields.params().max_solver_iters);
  EXPECT_LT(fields.irDropIterations(), fields.params().max_solver_iters);
  EXPECT_GT(fields.summary().max_temperature_c, fields.params().ambient_c);
  EXPECT_GT(fields.summary().max_ir_drop_v, 0.0f);
  EXPECT_FLOAT_EQ(fields.summary().max_stress, 1.0f);
  // Generous bound so a debug build still passes; Release is ~10x faster.
  EXPECT_LT(secs, 10.0);
}

}  // namespace
}  // namespace gpl
