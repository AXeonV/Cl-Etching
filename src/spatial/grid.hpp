#pragma once
#include "common/vector.hpp"
#include "common/material.hpp"
#include <vector>
#include <set>
namespace etch::spatial {
struct Geometry {
    int nx,ny,nz,surface_z,mask_thickness,opening,opening_y;
    int atoms_per_cell=30,deletion_threshold=5;
    double cell_nm=2.5;
};
struct Cell { int x,y,z; };
class Grid {
public:
    Geometry geometry;
    std::vector<uint8_t> material,boundary;
    std::vector<Sites> sites;
    std::set<int> surface_cells;
    explicit Grid(Geometry geometry);
    int wrap_x(int x) const { return (x%geometry.nx+geometry.nx)%geometry.nx; }
    int wrap_y(int y) const { return (y%geometry.ny+geometry.ny)%geometry.ny; }
    int index(int x,int y,int z) const { return (z*geometry.ny+wrap_y(y))*geometry.nx+wrap_x(x); }
    Cell coordinates(int cell) const { return {cell%geometry.nx,(cell/geometry.nx)%geometry.ny,cell/(geometry.nx*geometry.ny)}; }
    Material at(int x,int y,int z) const;
    int remaining(int cell) const { return site_count(sites[cell]); }
    void wrap_position(Vec3& p) const;
    void refresh(int x,int y,int z,int radius);
    void remove_empty_cell(int cell);
    long long silicon_inventory() const { return etch::silicon_inventory(sites); }
    long long chlorine_inventory() const { return etch::chlorine_inventory(sites); }
};
}
