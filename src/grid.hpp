#pragma once
#include "vector.hpp"
#include <array>
#include <cstdint>
#include <set>
#include <vector>

namespace etch {
enum Material : uint8_t { Gas, Silicon, Mask };
struct Geometry {
    int nx, nz;
    int surface_z, mask_thickness, opening;
    int atoms_per_cell = 30;
    int deletion_threshold = 5;
    double cell_nm = 2.5;
};
class Grid {
public:
    Geometry geometry;
    std::vector<uint8_t> material, boundary;
    std::set<int> surface_cells;
    std::vector<std::array<uint16_t, 4>> sites;
    explicit Grid(Geometry geometry);
    int wrap(int x) const { return (x % geometry.nx + geometry.nx) % geometry.nx; }
    int index(int x, int z) const { return z * geometry.nx + wrap(x); }
    Material at(int x, int z) const;
    int remaining(int index) const;
    void remove_empty_cell(int x, int z);
    void refresh(int x0, int z0, int radius);
    long long silicon_inventory() const;
    long long chlorine_inventory() const;
};
}

