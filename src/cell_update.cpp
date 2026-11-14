#include "cell_update.hpp"
#include <algorithm>

namespace etch {
CellUpdater::CellUpdater(Grid& g, Surface& s, double weight)
    : grid(g), surface(s), normal_weight(weight) {}

CellRemoval CellUpdater::remove(int x, int z) {
    int original = grid.index(x, z), selected = original;
    Curve curve = surface.curve(x, z);
    auto center = curve.project(curve.origin);
    Vec normal = curve.normal(center.u);
    // Solid-side distance is positive; gas-side protrusions have negative distance.
    double best = -center.distance;
    int selected_x = x, selected_z = z;
    for (int dz = -1; dz <= 1; ++dz) for (int dx = -1; dx <= 1; ++dx) {
        int zz = z + dz;
        if ((!dx && !dz) || zz < 0 || zz >= grid.geometry.nz) continue;
        int cell = grid.index(x + dx, zz);
        if (grid.material[cell] != Silicon || !grid.boundary[cell]) continue;
        Vec offset{static_cast<double>(dx), static_cast<double>(dz)};
        double angle = std::acos(std::clamp(dot(unit(offset), normal), -1.0, 1.0));
        double distance = curve.project(curve.origin + offset).distance;
        double score = normal_weight * angle - distance;
        if (score < best) {
            best = score; selected = cell;
            selected_x = grid.wrap(x + dx); selected_z = zz;
        }
    }
    auto residual = grid.sites[original];
    if (selected != original) grid.sites[original] = grid.sites[selected];
    grid.sites[selected] = {};
    grid.remove_empty_cell(selected_x, selected_z);
    surface.invalidate(selected_x, selected_z);
    return {selected_x, selected_z, residual, selected != original};
}
}
