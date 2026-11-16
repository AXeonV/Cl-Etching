#pragma once
#include "common/material.hpp"
#include "surface.hpp"

namespace etch::planar {
// Liu, section 3.6.2: select a removal location and transfer its inventory.
class CellUpdater {
    Grid& grid;
    Surface& surface;
    double normal_weight;
public:
    CellUpdater(Grid& grid, Surface& surface, double normal_weight = 0.3);
    CellRemoval remove(int cell);
};
}
