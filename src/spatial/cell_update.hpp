#pragma once
#include "common/cell_update.hpp"
#include "surface.hpp"

namespace etch::spatial {
class CellUpdater {
public:
  CellUpdateDiagnostics diagnostics;

private:
  Grid& grid;
  Surface& surface;
  double normal_weight;
  CellUpdate<Grid> update;
  void erase(int cell);

public:
  CellUpdater(
    Grid& grid, Surface& surface, double normal_weight = 0.3, unsigned long long seed = 42);
  CellRemoval remove(int cell, Vec3 collision_normal, Vec3 collision_position);
};
}
