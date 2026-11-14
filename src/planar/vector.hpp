#pragma once
#include <cmath>

namespace etch::planar {
struct Vec {
    double x = 0, z = 0;
    Vec operator+(Vec b) const { return {x + b.x, z + b.z}; }
    Vec operator-(Vec b) const { return {x - b.x, z - b.z}; }
    Vec operator*(double a) const { return {a * x, a * z}; }
};
inline double dot(Vec a, Vec b) { return a.x * b.x + a.z * b.z; }
inline double norm(Vec a) { return std::sqrt(dot(a, a)); }
inline Vec unit(Vec a) { return a * (1.0 / norm(a)); }
inline Vec reflected(Vec direction, Vec normal) {
    return direction - normal * (2 * dot(direction, normal));
}
}

