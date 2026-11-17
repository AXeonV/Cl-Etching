#include "backend.hpp"

namespace etch::planar {
void Backend::write_surface(Surface& surface, const Grid& grid, const std::string& path) {
  std::ofstream file(path);
  file << "x,z,tx,tz,nx,nz,a,b,c,cell\n" << std::setprecision(15);
  for (int cell : grid.surface_cells) {
    const auto& c = surface.curve(cell % grid.geometry.nx, cell / grid.geometry.nx);
    if (!c.defined)
      continue;
    file << c.origin.x << ',' << c.origin.z << ',' << c.tangent.x << ',' << c.tangent.z << ','
         << c.outward.x << ',' << c.outward.z << ',' << c.a << ',' << c.b << ',' << c.c << ','
         << cell << '\n';
  }
}

void Backend::write_field(Surface& surface, const Grid& grid, const std::string& path) {
  std::ofstream file(path, std::ios::binary);
  int n = surface.resolution();
  for (int z = 0; z <= grid.geometry.nz * n; ++z)
    for (int x = 0; x <= grid.geometry.nx * n; ++x) {
      auto node = surface.node(x, z);
      file.write(reinterpret_cast<const char*>(node.data()), sizeof(node));
    }
}

void Backend::write_geometry(Summary& out, const Geometry& g) {
  out.field("dimension", 2);
  out.field("nx", g.nx);
  out.field("nz", g.nz);
  out.field("cell_nm", g.cell_nm);
  out.field("surface_z", g.surface_z);
  out.field("opening_cells", g.opening);
  out.field("mask_cells", g.mask_thickness);
}
}
