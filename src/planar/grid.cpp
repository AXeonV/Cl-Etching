#include "grid.hpp"
#include <algorithm>

namespace etch::planar {
Grid::Grid(Geometry g)
    : geometry(g), material(g.nx * g.nz), boundary(g.nx * g.nz), sites(g.nx * g.nz) {
  int left = (g.nx - g.opening) / 2, right = left + g.opening;
  for (int z = 0; z < g.nz; ++z)
    for (int x = 0; x < g.nx; ++x) {
      int i = index(x, z);
      if (z >= g.surface_z) {
        material[i] = Silicon;
        sites[i][0] = g.atoms_per_cell;
      } else if (z >= g.surface_z - g.mask_thickness && (x < left || x >= right)) {
        material[i] = Mask;
      }
    }
  refresh(g.nx / 2, g.nz / 2, std::max(g.nx, g.nz));
}

std::array<int, 4> Grid::face_neighbors(int cell) const {
  int x = cell % geometry.nx, z = cell / geometry.nx;
  return {index(x - 1, z), index(x + 1, z), z > 0 ? index(x, z - 1) : -1,
    z + 1 < geometry.nz ? index(x, z + 1) : -1};
}

Material Grid::at(int x, int z) const {
  if (z < 0)
    return Gas;
  if (z >= geometry.nz)
    return Silicon;
  return static_cast<Material>(material[index(x, z)]);
}

void Grid::refresh(int x0, int z0, int radius) {
  int zmin = std::max(0, z0 - radius), zmax = std::min(geometry.nz - 1, z0 + radius);
  int xmin = x0 - std::min(radius, geometry.nx / 2);
  int xmax = radius >= geometry.nx / 2 ? xmin + geometry.nx - 1 : x0 + radius;
  for (int z = zmin; z <= zmax; ++z)
    for (int x = xmin; x <= xmax; ++x) {
      int i = index(x, z);
      boundary[i] = at(x, z) != Gas && (at(x - 1, z) == Gas || at(x + 1, z) == Gas ||
                                         at(x, z - 1) == Gas || at(x, z + 1) == Gas);
      if (boundary[i])
        surface_cells.insert(i);
      else
        surface_cells.erase(i);
    }
}

void Grid::remove_empty_cell(int x, int z) {
  material[index(x, z)] = Gas;
  refresh(x, z, 2);
}
}
