// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026-2026, The OpenROAD Authors

#include "physicsHeatMap.h"

#include <algorithm>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "gpl/Replace.h"
#include "odb/db.h"
#include "odb/geom.h"
#include "physicsBase.h"
#include "utl/Logger.h"
#include "web/core.h"
#include "web/heatMap.h"

namespace gpl {

namespace {

constexpr double kDefaultGridSizeUm = 10.0;

const std::vector<PhysicsHeatMapField>& allFields()
{
  static const std::vector<PhysicsHeatMapField> fields{
      PhysicsHeatMapField::kTemperature,
      PhysicsHeatMapField::kIrDrop,
      PhysicsHeatMapField::kStress,
      PhysicsHeatMapField::kPowerDensity};
  return fields;
}

std::vector<web::HeatMapSourceHandle>& registrations()
{
  static std::vector<web::HeatMapSourceHandle> handles;
  return handles;
}

// Instances handed to the Gui.  gpl registers its sources after Gui::init()
// has already instantiated the registered sources, so the instances are
// created here and kept alive for the lifetime of the process.
std::vector<std::unique_ptr<web::HeatMapDataSource>>& ownedInstances()
{
  static std::vector<std::unique_ptr<web::HeatMapDataSource>> instances;
  return instances;
}

}  // namespace

PhysicsHeatMapDataSource::PhysicsHeatMapDataSource(PhysicsHeatMapField field,
                                                   SnapshotProvider provider,
                                                   utl::Logger* logger)
    : web::RealValueHeatMapDataSource(logger,
                                      getFieldUnits(field),
                                      getFieldName(field),
                                      getFieldShortName(field),
                                      getFieldShortName(field)),
      field_(field),
      provider_(std::move(provider))
{
}

std::string PhysicsHeatMapDataSource::getFieldName(PhysicsHeatMapField field)
{
  switch (field) {
    case PhysicsHeatMapField::kTemperature:
      return "Physics: Temperature";
    case PhysicsHeatMapField::kIrDrop:
      return "Physics: IR Drop";
    case PhysicsHeatMapField::kStress:
      return "Physics: Reliability Stress";
    case PhysicsHeatMapField::kPowerDensity:
      return "Physics: Power Density";
  }
  return "Physics";
}

std::string PhysicsHeatMapDataSource::getFieldShortName(
    PhysicsHeatMapField field)
{
  switch (field) {
    case PhysicsHeatMapField::kTemperature:
      return "PhysicsTemperature";
    case PhysicsHeatMapField::kIrDrop:
      return "PhysicsIRDrop";
    case PhysicsHeatMapField::kStress:
      return "PhysicsStress";
    case PhysicsHeatMapField::kPowerDensity:
      return "PhysicsPowerDensity";
  }
  return "Physics";
}

std::string PhysicsHeatMapDataSource::getFieldUnits(PhysicsHeatMapField field)
{
  switch (field) {
    case PhysicsHeatMapField::kTemperature:
      return "C";
    case PhysicsHeatMapField::kIrDrop:
      return "V";
    case PhysicsHeatMapField::kStress:
      return "";
    case PhysicsHeatMapField::kPowerDensity:
      return "W/um^2";
  }
  return "";
}

const std::vector<float>* PhysicsHeatMapDataSource::getFieldValues(
    const PhysicsSnapshot& snap,
    PhysicsHeatMapField field)
{
  const std::vector<float>* values = nullptr;
  switch (field) {
    case PhysicsHeatMapField::kTemperature:
      values = &snap.temperature;
      break;
    case PhysicsHeatMapField::kIrDrop:
      values = &snap.ir_drop;
      break;
    case PhysicsHeatMapField::kStress:
      values = &snap.stress;
      break;
    case PhysicsHeatMapField::kPowerDensity:
      values = &snap.power_density;
      break;
  }
  if (values == nullptr || values->empty()) {
    return nullptr;
  }
  return values;
}

const PhysicsSnapshot* PhysicsHeatMapDataSource::getSnapshot() const
{
  if (!provider_) {
    return nullptr;
  }
  const PhysicsSnapshot* snap = provider_();
  if (snap == nullptr || snap->grid.cnt_x <= 0 || snap->grid.cnt_y <= 0
      || snap->bin_w_dbu <= 0 || snap->bin_h_dbu <= 0) {
    return nullptr;
  }
  return snap;
}

double PhysicsHeatMapDataSource::getDbuPerMicron() const
{
  odb::dbBlock* block = getBlock();
  if (block == nullptr) {
    return 1.0;
  }
  return block->getDbUnitsPerMicron();
}

double PhysicsHeatMapDataSource::getGridXSize() const
{
  const PhysicsSnapshot* snap = getSnapshot();
  if (snap == nullptr) {
    return kDefaultGridSizeUm;
  }
  return snap->bin_w_dbu / getDbuPerMicron();
}

double PhysicsHeatMapDataSource::getGridYSize() const
{
  const PhysicsSnapshot* snap = getSnapshot();
  if (snap == nullptr) {
    return kDefaultGridSizeUm;
  }
  return snap->bin_h_dbu / getDbuPerMicron();
}

void PhysicsHeatMapDataSource::populateXYGrid()
{
  const PhysicsSnapshot* snap = getSnapshot();
  if (snap == nullptr) {
    web::RealValueHeatMapDataSource::populateXYGrid();
    return;
  }

  // Use the snapshot bins directly so every map cell is exactly one physics
  // bin and no resampling is needed.
  std::vector<int> x_grid;
  x_grid.reserve(snap->grid.cnt_x + 1);
  for (int x = 0; x <= snap->grid.cnt_x; x++) {
    x_grid.push_back(snap->grid_lx_dbu + x * snap->bin_w_dbu);
  }
  std::vector<int> y_grid;
  y_grid.reserve(snap->grid.cnt_y + 1);
  for (int y = 0; y <= snap->grid.cnt_y; y++) {
    y_grid.push_back(snap->grid_ly_dbu + y * snap->bin_h_dbu);
  }
  setXYMapGrid(x_grid, y_grid);
}

bool PhysicsHeatMapDataSource::populateMap()
{
  const PhysicsSnapshot* snap = getSnapshot();
  if (snap == nullptr) {
    return false;
  }
  const std::vector<float>* values = getFieldValues(*snap, field_);
  if (values == nullptr) {
    return false;
  }

  const int cnt_x = snap->grid.cnt_x;
  const int cnt_y = snap->grid.cnt_y;
  if (values->size() < static_cast<size_t>(cnt_x) * cnt_y) {
    return false;
  }

  for (int x = 0; x < cnt_x; x++) {
    const int lx = snap->grid_lx_dbu + x * snap->bin_w_dbu;
    const int ux = lx + snap->bin_w_dbu;
    for (int y = 0; y < cnt_y; y++) {
      const int ly = snap->grid_ly_dbu + y * snap->bin_h_dbu;
      const int uy = ly + snap->bin_h_dbu;
      const odb::Rect bin(lx, ly, ux, uy);
      addToMap(bin, (*values)[snap->grid.index(x, y)]);
    }
  }
  return true;
}

void PhysicsHeatMapDataSource::combineMapData(bool base_has_value,
                                              double& base,
                                              const double new_data,
                                              const double data_area,
                                              const double intersection_area,
                                              const double rect_area)
{
  // Bins map one-to-one onto the grid cells; ignore neighbours that merely
  // touch a cell and keep the worst (largest) value otherwise.
  if (intersection_area <= 0.0) {
    return;
  }
  if (!base_has_value) {
    base = new_data;
  } else {
    base = std::max(base, new_data);
  }
}

void registerPhysicsHeatMaps(Replace* replace, utl::Logger* logger)
{
  if (!registrations().empty()) {
    return;
  }
  for (const PhysicsHeatMapField field : allFields()) {
    auto factory = [replace, logger, field]() {
      return std::make_unique<PhysicsHeatMapDataSource>(
          field,
          [replace]() -> const PhysicsSnapshot* {
            return replace ? replace->getPhysicsSnapshot() : nullptr;
          },
          logger);
    };
    registrations().push_back(web::registerHeatMapSource(
        PhysicsHeatMapDataSource::getFieldName(field),
        PhysicsHeatMapDataSource::getFieldShortName(field),
        PhysicsHeatMapDataSource::getFieldShortName(field),
        factory));
    auto instance = factory();
    web::Gui::get()->registerHeatMap(instance.get());
    ownedInstances().push_back(std::move(instance));
  }
}

void invalidatePhysicsHeatMaps()
{
  for (auto& handle : registrations()) {
    if (handle) {
      handle->invalidateInstances();
    }
  }
}

}  // namespace gpl
