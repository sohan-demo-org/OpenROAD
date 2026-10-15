// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026-2026, The OpenROAD Authors

#pragma once

namespace utl {
class Logger;
}

namespace gpl {

struct PhysicsSnapshot;

// Prints the `report_physics_placement` summary and per-pass history table
// from the snapshot left by `global_placement -physics_driven`.  Warns
// (GPL-0230) and returns when there is no physics data.
void reportPhysicsPlacement(const PhysicsSnapshot* snapshot,
                            utl::Logger* logger);

}  // namespace gpl
