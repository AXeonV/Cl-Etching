#pragma once
#include "surface.hpp"
namespace etch::spatial {
struct CellRemoval { Sites residual;bool relocated; };
class CellUpdater {
    Grid& grid;
    Surface& surface;
    double normal_weight;
public:
    CellUpdater(Grid& grid,Surface& surface,double normal_weight=0.3):grid(grid),surface(surface),normal_weight(normal_weight) {}
    CellRemoval remove(int cell);
};
}
