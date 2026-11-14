#include "surface.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

namespace etch {
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
        if (t >= start && t <= end && dot(d, normal(u + du * t)) < -1e-9)
            nearest = std::min(nearest, t);
    }
    if (std::isfinite(nearest)) return nearest;
    return std::nullopt;
}

// Symmetric eigendecomposition gives the Moore-Penrose least-squares solution.
// It also exposes the rank of the local fitting matrix.
int least_squares(double matrix[3][3], const double rhs[3], double result[3]) {
    double eigenvectors[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    for (int iteration = 0; iteration < 40; ++iteration) {
        int p = 0, q = 1;
        for (int i = 0; i < 3; ++i) for (int j = i + 1; j < 3; ++j)
            if (std::abs(matrix[i][j]) > std::abs(matrix[p][q])) { p = i; q = j; }
        double scale = std::abs(matrix[0][0]) + std::abs(matrix[1][1]) + std::abs(matrix[2][2]);
        if (std::abs(matrix[p][q]) < 1e-14 * scale) break;
        double angle = 0.5 * std::atan2(2 * matrix[p][q], matrix[q][q] - matrix[p][p]);
        double c = std::cos(angle), s = std::sin(angle);
        double pp = matrix[p][p], qq = matrix[q][q], pq = matrix[p][q];
        for (int k = 0; k < 3; ++k) if (k != p && k != q) {
            double kp = matrix[k][p], kq = matrix[k][q];
            matrix[k][p] = matrix[p][k] = c * kp - s * kq;
            matrix[k][q] = matrix[q][k] = s * kp + c * kq;
        }
        matrix[p][p] = c*c*pp - 2*c*s*pq + s*s*qq;
        matrix[q][q] = s*s*pp + 2*c*s*pq + c*c*qq;
        matrix[p][q] = matrix[q][p] = 0;
        for (int k = 0; k < 3; ++k) {
            double vp = eigenvectors[k][p], vq = eigenvectors[k][q];
            eigenvectors[k][p] = c * vp - s * vq;
            eigenvectors[k][q] = s * vp + c * vq;
        }
    }
    double maximum = std::max({matrix[0][0], matrix[1][1], matrix[2][2]});
    int rank = 0;
    for (int j = 0; j < 3; ++j) if (matrix[j][j] > maximum * 1e-12) {
        ++rank;
        double projected = 0;
        for (int k = 0; k < 3; ++k) projected += eigenvectors[k][j] * rhs[k];
        for (int k = 0; k < 3; ++k) result[k] += eigenvectors[k][j] * projected / matrix[j][j];
    }
    return rank;
}

Surface::Surface(Grid& g, int r) : grid(g), radius(r), curves(g.material.size()), dirty(g.material.size(), 1) {}
const Curve& Surface::curve(int x, int z) {
    int id = grid.index(x, z);
    if (!dirty[id]) return curves[id];
    std::vector<Vec> points;
    Vec gas_direction{};
    for (int dz = -radius; dz <= radius; ++dz) for (int dx = -radius; dx <= radius; ++dx) {
        int zz = z + dz;
        if (zz < 0 || zz >= grid.geometry.nz || !grid.boundary[grid.index(x + dx, zz)]) continue;
        points.push_back({static_cast<double>(dx), static_cast<double>(dz)});
        for (Vec n : {Vec{-1, 0}, Vec{1, 0}, Vec{0, -1}, Vec{0, 1}})
            if (grid.at(x + dx + static_cast<int>(n.x), zz + static_cast<int>(n.z)) == Gas)
                gas_direction = gas_direction + n;
    }
    Vec mean{};
    for (Vec p : points) mean = mean + p;
    mean = mean * (1.0 / points.size());
    double xx = 0, xz = 0, zz = 0;
    for (Vec p : points) { p = p - mean; xx += p.x*p.x; xz += p.x*p.z; zz += p.z*p.z; }
    double angle = 0.5 * std::atan2(2 * xz, xx - zz);
    Vec tangent{std::cos(angle), std::sin(angle)}, outward{tangent.z, -tangent.x};
    if (dot(outward, gas_direction) < 0) outward = outward * -1;
    double matrix[3][3] = {}, rhs[3] = {}, coefficients[3] = {};
    for (Vec p : points) {
        double u = dot(p, tangent), v = dot(p, outward), row[] = {u*u, u, 1};
        for (int j = 0; j < 3; ++j) {
            rhs[j] += row[j] * v;
            for (int k = 0; k < 3; ++k) matrix[j][k] += row[j] * row[k];
        }
    }
    int rank = least_squares(matrix, rhs, coefficients);
    if (rank < 3) ++rank_deficient_fits;
    curves[id] = {{x + 0.5, z + 0.5}, tangent, outward, coefficients[0], coefficients[1], coefficients[2]};
    dirty[id] = 0;
    ++fit_count;
    return curves[id];
}
void Surface::invalidate(int x, int z) {
    for (int dz = -radius - 1; dz <= radius + 1; ++dz)
        for (int dx = -radius - 1; dx <= radius + 1; ++dx)
            if (z + dz >= 0 && z + dz < grid.geometry.nz) dirty[grid.index(x + dx, z + dz)] = 1;
}
}
