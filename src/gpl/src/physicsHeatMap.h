// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026-2026, The OpenROAD Authors

#pragma once

#include <functional>
#include <string>
#include <vector>

#include "web/heatMap.h"

namespace utl {
class Logger;
}

namespace gpl {

class Replace;
struct PhysicsSnapshot;

// Which field of the PhysicsSnapshot a heat map renders.
enum class PhysicsHeatMapField
{
  kTemperature,
  kIrDrop,
  kStress,
  kPowerDensity
};

// GUI heat map backed by the PhysicsSnapshot left behind by
// `global_placement -physics_driven`.  The snapshot is fetched lazily through
// the provider every time the map is (re)built, so the data source is valid
// before any placement has run (it then renders nothing) and keeps working
// after `gpl::replace_reset_cmd`.
class PhysicsHeatMapDataSource : public web::RealValueHeatMapDataSource
{
 public:
  using SnapshotProvider = std::function<const PhysicsSnapshot*()>;

  PhysicsHeatMapDataSource(PhysicsHeatMapField field,
                           SnapshotProvider provider,
                           utl::Logger* logger);

  static std::string getFieldName(PhysicsHeatMapField field);
  static std::string getFieldShortName(PhysicsHeatMapField field);
  static std::string getFieldUnits(PhysicsHeatMapField field);
  // Returns nullptr when the snapshot does not carry this field.
  static const std::vector<float>* getFieldValues(const PhysicsSnapshot& snap,
                                                  PhysicsHeatMapField field);

  PhysicsHeatMapField getField() const { return field_; }

  bool canAdjustGrid() const override { return false; }
  double getGridXSize() const override;
  double getGridYSize() const override;

 protected:
  void populateXYGrid() override;
  bool populateMap() override;
  void combineMapData(bool base_has_value,
                      double& base,
                      double new_data,
                      double data_area,
                      double intersection_area,
                      double rect_area) override;

 private:
  // Snapshot with a usable grid, or nullptr.
  const PhysicsSnapshot* getSnapshot() const;
  double getDbuPerMicron() const;

  PhysicsHeatMapField field_;
  SnapshotProvider provider_;
};

// Registers the four physics heat maps ("Physics: Temperature", "Physics: IR
// Drop", "Physics: Reliability Stress", "Physics: Power Density") with the GUI.
// Safe to call without a GUI; the sources are instantiated when one starts.
void registerPhysicsHeatMaps(Replace* replace, utl::Logger* logger);

// Drops the cached maps so the next paint re-reads the snapshot.  Call after a
// placement run has produced a new snapshot.
void invalidatePhysicsHeatMaps();

}  // namespace gpl
