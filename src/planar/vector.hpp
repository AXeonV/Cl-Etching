#pragma once
#include "common/vector.hpp"

namespace etch::planar {
struct Vec {
    double x = 0, z = 0;
    double operator[](int axis) const { return axis==0?x:z; }
    Vec operator+(Vec b) const { return {x + b.x, z + b.z}; }
    Vec operator-(Vec b) const { return {x - b.x, z - b.z}; }
    Vec operator*(double a) const { return {a * x, a * z}; }
};
inline double dot(Vec a, Vec b) { return a.x * b.x + a.z * b.z; }
inline double norm(Vec a) { return std::sqrt(dot(a, a)); }
inline Vec unit(Vec a) { return a * (1.0 / norm(a)); }
inline Vec projection(Vec3 direction) { return {direction.x, direction.z}; }
}
