#pragma once
#include "surface.hpp"

namespace etch::planar {
struct CellRemoval {
    int x, z;
    std::array<uint16_t, 4> residual;
    bool relocated;
};

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
