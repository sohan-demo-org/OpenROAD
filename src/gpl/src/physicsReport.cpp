// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026-2026, The OpenROAD Authors

#include "physicsReport.h"

#include <algorithm>
#include <cstddef>
#include <vector>

#include "physicsBase.h"
#include "utl/Logger.h"

namespace gpl {

namespace {

bool hasPhysicsData(const PhysicsSnapshot* snap)
{
  return snap != nullptr && snap->grid.cnt_x > 0 && snap->grid.cnt_y > 0
         && !snap->temperature.empty();
}

// Flat index of the hottest bin, or -1 when the field is empty.
int hottestBin(const PhysicsSnapshot& snap)
{
  if (snap.temperature.empty()) {
    return -1;
  }
  const auto it
      = std::max_element(snap.temperature.begin(), snap.temperature.end());
  return static_cast<int>(it - snap.temperature.begin());
}

}  // namespace

void reportPhysicsPlacement(const PhysicsSnapshot* snap, utl::Logger* logger)
{
  if (!hasPhysicsData(snap)) {
    logger->warn(utl::GPL,
                 230,
                 "no physics data; run global_placement -physics_driven first");
    return;
  }

  const PhysicsSummary& s = snap->summary;

  logger->report("");
  logger->report("Physics-driven placement report");
  logger->report("--------------------------------------------------");
  logger->report("Power source           : {}",
                 snap->used_sta_power ? "OpenSTA (liberty power)"
                                      : "proxy (leakage / cell area)");
  logger->report("Total power            : {:.4e} W", s.total_power_w);
  logger->report("Max temperature        : {:.2f} C", s.max_temperature_c);
  logger->report("Avg temperature        : {:.2f} C", s.avg_temperature_c);
  logger->report("Max IR drop            : {:.3f} mV ({:.2f}%)",
                 1000.0 * s.max_ir_drop_v,
                 s.max_ir_drop_percent);
  logger->report("Max delay derate       : {:.4f}", s.max_delay_derate);
  logger->report("Max reliability stress : {:.4f}", s.max_stress);
  logger->report("Physics passes         : {}", snap->passes);
  logger->report("Inflated area          : {:.2f}%",
                 snap->inflated_area_percent);

  const int hot = hottestBin(*snap);
  if (hot >= 0) {
    const int hx = hot / snap->grid.cnt_y;
    const int hy = hot % snap->grid.cnt_y;
    const double x_um = snap->grid.lx_um + (hx + 0.5) * snap->grid.bin_w_um;
    const double y_um = snap->grid.ly_um + (hy + 0.5) * snap->grid.bin_h_um;
    logger->report("Hottest bin            : ({}, {}) at ({:.2f}, {:.2f}) um",
                   hx,
                   hy,
                   x_um,
                   y_um);
  }
  logger->report("Grid                   : {} x {} bins of {:.2f} x {:.2f} um",
                 snap->grid.cnt_x,
                 snap->grid.cnt_y,
                 snap->grid.bin_w_um,
                 snap->grid.bin_h_um);

  if (!snap->history.empty()) {
    logger->report("");
    logger->report("Pass history");
    logger->report("{:>4} {:>6} {:>9} {:>9} {:>10} {:>10}",
                   "pass",
                   "iter",
                   "overflow",
                   "Tmax(C)",
                   "IRmax(mV)",
                   "derateMax");
    int pass = 1;
    for (const auto& rec : snap->history) {
      logger->report("{:>4} {:>6} {:>9.4f} {:>9.2f} {:>10.3f} {:>10.4f}",
                     pass++,
                     rec.iter,
                     rec.overflow,
                     rec.summary.max_temperature_c,
                     1000.0 * rec.summary.max_ir_drop_v,
                     rec.summary.max_delay_derate);
    }
  }
  logger->report("");
}

}  // namespace gpl
