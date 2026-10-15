// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026-2026, The OpenROAD Authors

#include <algorithm>
#include <cstddef>
#include <memory>
#include <vector>

#include "gtest/gtest.h"
#include "odb/db.h"
#include "odb/geom.h"
#include "src/gpl/src/physicsBase.h"
#include "src/gpl/src/physicsHeatMap.h"
#include "tst/db_fixture.h"

namespace gpl {
namespace {

constexpr int kDbuPerMicron = 1000;
constexpr int kCntX = 3;
constexpr int kCntY = 2;
constexpr int kBinDbu = 10 * kDbuPerMicron;

class PhysicsHeatMapTest : public tst::DbFixture
{
 protected:
  PhysicsHeatMapTest()
  {
    db_->setDbuPerMicron(kDbuPerMicron);
    tech_ = odb::dbTech::create(db_.get(), "tech");
    chip_ = odb::dbChip::create(db_.get(), tech_, "chip");
    block_ = odb::dbBlock::create(chip_, "top");
    block_->setDieArea(odb::Rect(0, 0, kCntX * kBinDbu, kCntY * kBinDbu));
  }

  // 3 x 2 grid where every field value encodes its (x, y) bin as x * 10 + y.
  PhysicsSnapshot makeSnapshot() const
  {
    PhysicsSnapshot snap;
    snap.grid.cnt_x = kCntX;
    snap.grid.cnt_y = kCntY;
    snap.grid.bin_w_um = 10;
    snap.grid.bin_h_um = 10;
    snap.bin_w_dbu = kBinDbu;
    snap.bin_h_dbu = kBinDbu;
    const size_t n = static_cast<size_t>(kCntX) * kCntY;
    snap.temperature.resize(n);
    snap.ir_drop.resize(n);
    snap.stress.resize(n);
    snap.power_density.resize(n);
    for (int x = 0; x < kCntX; x++) {
      for (int y = 0; y < kCntY; y++) {
        const float v = x * 10 + y;
        const int i = snap.grid.index(x, y);
        snap.temperature[i] = 25 + v;
        snap.ir_drop[i] = 1e-3f * v;
        snap.stress[i] = 0.01f * v;
        snap.power_density[i] = 1e-6f * v;
      }
    }
    return snap;
  }

  std::unique_ptr<PhysicsHeatMapDataSource> makeSource(
      PhysicsHeatMapField field,
      const PhysicsSnapshot* snap)
  {
    auto source = std::make_unique<PhysicsHeatMapDataSource>(
        field, [snap]() { return snap; }, getLogger());
    source->setChip(chip_);
    return source;
  }

  odb::dbTech* tech_ = nullptr;
  odb::dbChip* chip_ = nullptr;
  odb::dbBlock* block_ = nullptr;
};

TEST_F(PhysicsHeatMapTest, NamesAndUnits)
{
  EXPECT_EQ(
      PhysicsHeatMapDataSource::getFieldName(PhysicsHeatMapField::kTemperature),
      "Physics: Temperature");
  EXPECT_EQ(PhysicsHeatMapDataSource::getFieldUnits(
                PhysicsHeatMapField::kTemperature),
            "C");
  EXPECT_EQ(
      PhysicsHeatMapDataSource::getFieldName(PhysicsHeatMapField::kIrDrop),
      "Physics: IR Drop");
  EXPECT_EQ(
      PhysicsHeatMapDataSource::getFieldUnits(PhysicsHeatMapField::kIrDrop),
      "V");
  EXPECT_EQ(
      PhysicsHeatMapDataSource::getFieldName(PhysicsHeatMapField::kStress),
      "Physics: Reliability Stress");
  EXPECT_EQ(
      PhysicsHeatMapDataSource::getFieldUnits(PhysicsHeatMapField::kStress),
      "");
  EXPECT_EQ(PhysicsHeatMapDataSource::getFieldName(
                PhysicsHeatMapField::kPowerDensity),
            "Physics: Power Density");
  EXPECT_EQ(PhysicsHeatMapDataSource::getFieldUnits(
                PhysicsHeatMapField::kPowerDensity),
            "W/um^2");
}

TEST_F(PhysicsHeatMapTest, EmptySnapshotProducesNoMap)
{
  // No snapshot at all (before any -physics_driven run).
  {
    auto source = makeSource(PhysicsHeatMapField::kTemperature, nullptr);
    EXPECT_FALSE(source->canAdjustGrid());
    source->ensureMap();
    EXPECT_FALSE(source->isPopulated());
  }
  // Snapshot exists but its grid / fields are empty (stubbed coupling).
  {
    PhysicsSnapshot empty;
    auto source = makeSource(PhysicsHeatMapField::kIrDrop, &empty);
    source->ensureMap();
    EXPECT_FALSE(source->isPopulated());
  }
  // Grid is set but the requested field is missing.
  {
    PhysicsSnapshot snap = makeSnapshot();
    snap.stress.clear();
    auto source = makeSource(PhysicsHeatMapField::kStress, &snap);
    source->ensureMap();
    EXPECT_FALSE(source->isPopulated());
  }
}

TEST_F(PhysicsHeatMapTest, GridSizesFollowSnapshotBins)
{
  const PhysicsSnapshot snap = makeSnapshot();
  auto source = makeSource(PhysicsHeatMapField::kTemperature, &snap);
  EXPECT_DOUBLE_EQ(source->getGridXSize(), 10.0);
  EXPECT_DOUBLE_EQ(source->getGridYSize(), 10.0);
}

TEST_F(PhysicsHeatMapTest, MapMatchesSnapshotBins)
{
  const PhysicsSnapshot snap = makeSnapshot();
  const struct
  {
    PhysicsHeatMapField field;
    const std::vector<float>* values;
  } cases[] = {{PhysicsHeatMapField::kTemperature, &snap.temperature},
               {PhysicsHeatMapField::kIrDrop, &snap.ir_drop},
               {PhysicsHeatMapField::kStress, &snap.stress},
               {PhysicsHeatMapField::kPowerDensity, &snap.power_density}};

  for (const auto& c : cases) {
    auto source = makeSource(c.field, &snap);
    source->ensureMap();
    ASSERT_TRUE(source->isPopulated());

    // The map holds each value as a percentage of the field's [min, max]
    // range (rounded to 1e-3).
    const auto [min_it, max_it]
        = std::minmax_element(c.values->begin(), c.values->end());
    const double range = *max_it - *min_it;

    const auto& map = source->getMap();
    ASSERT_EQ(map.shape()[0], static_cast<size_t>(kCntX));
    ASSERT_EQ(map.shape()[1], static_cast<size_t>(kCntY));
    for (int x = 0; x < kCntX; x++) {
      for (int y = 0; y < kCntY; y++) {
        const auto& cell = map[x][y];
        ASSERT_NE(cell, nullptr);
        EXPECT_TRUE(cell->has_value);
        const double expected_percent
            = 100.0 * ((*c.values)[snap.grid.index(x, y)] - *min_it) / range;
        EXPECT_NEAR(cell->value, expected_percent, 1e-2);
        EXPECT_EQ(cell->rect,
                  odb::Rect(x * kBinDbu,
                            y * kBinDbu,
                            (x + 1) * kBinDbu,
                            (y + 1) * kBinDbu));
      }
    }
  }
}

}  // namespace
}  // namespace gpl
