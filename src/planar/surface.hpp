#pragma once
#include "grid.hpp"
#include "common/implicit.hpp"
#include <unordered_map>

namespace etch::planar {
struct CurveProjection { double u, distance; };

// Local least-squares graph: v = a*u^2 + b*u + c in a rotated frame.
struct Curve {
    Vec origin, tangent, outward;
    double a = 0, b = 0, c = 0;
    bool defined = true;
    Vec position(double u) const;
    Vec normal(double u) const;
    CurveProjection project(Vec point) const;
};
class Surface {
    Grid& grid;
    int radius,subdivisions;
    using Field=ImplicitCell<2>;
    std::vector<Field::Node> node_values;
    std::unordered_map<int,Field> fields;
    std::vector<uint8_t> dirty_nodes;
    std::vector<Curve> curves;
    std::vector<uint8_t> dirty, dirty_candidates;
    std::vector<std::vector<int>> nearby_cells;
public:
    long long fit_count = 0, rank_deficient_fits = 0, undefined_normal_fits = 0;
    explicit Surface(Grid& grid, int radius = 3,int subdivisions = 1);
    Field::Node node(int x,int z);
    const Field& field(int x,int z);
    int resolution() const { return subdivisions; }
    int fit_radius() const { return radius; }
    const Curve& curve(int x, int z);
    const std::vector<int>& candidates(int x, int z);
    void invalidate(int x, int z);
};
}
