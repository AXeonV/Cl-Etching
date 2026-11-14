#include "surface.hpp"
#include "common/linear.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

namespace etch::planar {
Vec Curve::position(double u) const {
    return origin + tangent * u + outward * (a * u * u + b * u + c);
}
Vec Curve::normal(double u) const { return unit(outward - tangent * (2 * a * u + b)); }
// The nearest point solves (u-u0) + (f(u)-v0)*f'(u) = 0.
CurveProjection Curve::project(Vec point) const {
    Vec local = point - origin;
    double u0 = dot(local, tangent), v0 = dot(local, outward);
    double u = (u0 + b * (v0 - c)) / (1 + b*b);
    if (a != 0) {
        auto slope_equation = [&](double x) {
            return (x-u0) + (a*x*x+b*x+c-v0)*(2*a*x+b);
        };
        // A closer point cannot lie farther horizontally than this known distance.
        double bound = std::abs(a*u0*u0+b*u0+c-v0);
        std::vector<double> intervals{u0-bound, u0+bound};
        double A = 6*a*a, B = 6*a*b, C = 1+b*b+2*a*(c-v0);
        double discriminant = B*B - 4*A*C;
        if (discriminant > 0) {
            double q = -0.5*(B + std::copysign(std::sqrt(discriminant), B));
            for (double critical : {q/A, C/q})
                if (critical > intervals[0] && critical < intervals[1]) intervals.push_back(critical);
        }
        std::sort(intervals.begin(), intervals.end());
        u = u0;
        double nearest = norm(point - position(u));
        for (size_t i = 1; i < intervals.size(); ++i) {
            double left = intervals[i-1], right = intervals[i];
            double gl = slope_equation(left), gr = slope_equation(right);
            if ((gl > 0 && gr > 0) || (gl < 0 && gr < 0)) continue;
            if (gl == 0) right = left;
            else if (gr == 0) left = right;
            for (int step = 0; step < 60; ++step) {
                double middle = (left+right)/2, gm = slope_equation(middle);
                if (gm == 0) { left = right = middle; break; }
                if ((gl <= 0 && gm <= 0) || (gl >= 0 && gm >= 0)) { left = middle; gl = gm; }
                else right = middle;
            }
            double candidate = (left+right)/2, distance = norm(point-position(candidate));
            if (distance < nearest) { nearest = distance; u = candidate; }
        }
    }
    Vec delta = point-position(u);
    return {u, std::copysign(norm(delta), dot(delta, normal(u)))};
}

std::optional<double> Curve::intersection(Vec p, Vec d, double start, double end) const {
    Vec local = p - origin;
    double u = dot(local, tangent), v = dot(local, outward);
    double du = dot(d, tangent), dv = dot(d, outward);
    double A = a * du * du, B = 2 * a * u * du + b * du - dv;
    double C = a * u * u + b * u + c - v;
    double roots[2]; int count = 0;
    if (std::abs(A) < 1e-12) {
        if (std::abs(B) > 1e-12) roots[count++] = -C / B;
    } else {
        double discriminant = B * B - 4 * A * C;
        if (discriminant >= 0) {
            double q = -0.5 * (B + std::copysign(std::sqrt(discriminant), B));
            roots[count++] = q / A;
            if (std::abs(q) > 1e-15) roots[count++] = C / q;
        }
    }
    double nearest = std::numeric_limits<double>::infinity();
    for (int i = 0; i < count; ++i) {
        double t = roots[i];
        double hit_u = u + du*t;
        if (t >= start && t <= end && hit_u >= u_min && hit_u <= u_max
            && dot(d, normal(hit_u)) < -1e-9)
            nearest = std::min(nearest, t);
    }
    if (std::isfinite(nearest)) return nearest;
    return std::nullopt;
}

Surface::Surface(Grid& g, int r) : grid(g), radius(r), curves(g.material.size()), dirty(g.material.size(), 1),
    dirty_candidates(g.material.size(),1), nearby_cells(g.material.size()) {}
const Curve& Surface::curve(int x, int z) {
    int id = grid.index(x, z);
    if (!dirty[id]) return curves[id];
    std::vector<Vec> points;
    Vec solid_centroid{};
    int solid_count = 0;
    for (int dz = -radius; dz <= radius; ++dz) for (int dx = -radius; dx <= radius; ++dx) {
        int zz = z+dz;
        if (zz < 0 || zz >= grid.geometry.nz) continue;
        int cell = grid.index(x+dx,zz);
        if (grid.material[cell] != Gas) {
            solid_centroid = solid_centroid+Vec{static_cast<double>(dx),static_cast<double>(dz)};
            ++solid_count;
        }
        if (grid.boundary[cell]) points.push_back({static_cast<double>(dx),static_cast<double>(dz)});
    }
    solid_centroid = solid_centroid*(1.0/solid_count);
    if (solid_centroid.x==0 && solid_centroid.z==0) {
        curves[id] = {}; curves[id].origin = {x+0.5,z+0.5}; curves[id].defined = false;
        dirty[id] = 0; ++fit_count; ++undefined_normal_fits;
        return curves[id];
    }
    Vec outward = unit(solid_centroid*-1);
    Vec tangent{-outward.z,outward.x};
    double matrix[3][3] = {}, rhs[3] = {}, coefficients[3] = {};
    double u_min = std::numeric_limits<double>::infinity(), u_max = -u_min;
    for (Vec p : points) {
        double u = dot(p, tangent), v = dot(p, outward), row[] = {u*u, u, 1};
        u_min = std::min(u_min,u); u_max = std::max(u_max,u);
        for (int j = 0; j < 3; ++j) {
            rhs[j] += row[j] * v;
            for (int k = 0; k < 3; ++k) matrix[j][k] += row[j] * row[k];
        }
    }
    int rank = least_squares(matrix, rhs, coefficients);
    if (rank < 3) ++rank_deficient_fits;
    curves[id] = {{x+0.5,z+0.5},tangent,outward,
                  coefficients[0],coefficients[1],coefficients[2],u_min,u_max};
    dirty[id] = 0;
    ++fit_count;
    return curves[id];
}
const std::vector<int>& Surface::candidates(int x, int z) {
    int id = grid.index(x,z);
    if (!dirty_candidates[id]) return nearby_cells[id];
    auto& cells = nearby_cells[id]; cells.clear();
    int reach = radius+1;
    for (int dz=-reach;dz<=reach;++dz) for(int dx=-reach;dx<=reach;++dx) {
        int zz=z+dz;
        if (zz>=0 && zz<grid.geometry.nz && grid.boundary[grid.index(x+dx,zz)])
            cells.push_back(grid.index(x+dx,zz));
    }
    dirty_candidates[id]=0;
    return cells;
}

void Surface::invalidate(int x, int z) {
    for (int dz=-radius-2;dz<=radius+2;++dz) for(int dx=-radius-2;dx<=radius+2;++dx) {
        if (z+dz<0 || z+dz>=grid.geometry.nz) continue;
        int id=grid.index(x+dx,z+dz);
        dirty_candidates[id]=1;
        if (std::abs(dx)<=radius+1 && std::abs(dz)<=radius+1) dirty[id]=1;
    }
}
}
