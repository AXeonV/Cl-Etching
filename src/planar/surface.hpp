#pragma once
#include "grid.hpp"
#include <optional>
#include <limits>

namespace etch::planar {
struct CurveProjection { double u, distance; };

// Local least-squares graph: v = a*u^2 + b*u + c in a rotated frame.
struct Curve {
    Vec origin, tangent, outward;
    double a = 0, b = 0, c = 0;
    double u_min = -std::numeric_limits<double>::infinity();
    double u_max = std::numeric_limits<double>::infinity();
    bool defined = true;
    Vec position(double u) const;
    Vec normal(double u) const;
    CurveProjection project(Vec point) const;
    std::optional<double> intersection(Vec p, Vec d, double start, double end) const;
};
class Surface {
    Grid& grid;
    int radius;
    std::vector<Curve> curves;
    std::vector<uint8_t> dirty, dirty_candidates;
    std::vector<std::vector<int>> nearby_cells;
public:
    long long fit_count = 0, rank_deficient_fits = 0, undefined_normal_fits = 0;
    explicit Surface(Grid& grid, int radius = 3);
    int fit_radius() const { return radius; }
    const Curve& curve(int x, int z);
    const std::vector<int>& candidates(int x, int z);
    void invalidate(int x, int z);
};
}
