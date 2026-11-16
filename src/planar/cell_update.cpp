#include "cell_update.hpp"
#include <algorithm>

namespace etch::planar {
CellUpdater::CellUpdater(Grid& g, Surface& s, double weight)
    : grid(g), surface(s), normal_weight(weight) {}

CellRemoval CellUpdater::remove(int cell_id) {
    int x=cell_id%grid.geometry.nx, z=cell_id/grid.geometry.nx;
    int original = grid.index(x, z), selected = original;
    Curve curve = surface.curve(x, z);
    auto center = curve.project(curve.origin);
    Vec normal = curve.normal(center.u);
    // Solid-side distance is positive; gas-side protrusions have negative distance.
    double best = -center.distance;
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
        }
    }
    auto removal=transfer_sites_for_removal(grid.sites,original,selected);
    int removed_x=selected%grid.geometry.nx,removed_z=selected/grid.geometry.nx;
    grid.remove_empty_cell(removed_x,removed_z);
    surface.invalidate(removed_x,removed_z);
    return removal;
}
}
