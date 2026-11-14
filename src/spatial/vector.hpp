#pragma once
#include <cmath>
namespace etch::spatial {
struct Vec {
    double x=0,y=0,z=0;
    double operator[](int axis) const { return axis==0?x:axis==1?y:z; }
};
inline Vec operator+(Vec a,Vec b) { return {a.x+b.x,a.y+b.y,a.z+b.z}; }
inline Vec operator-(Vec a,Vec b) { return {a.x-b.x,a.y-b.y,a.z-b.z}; }
inline Vec operator*(Vec a,double s) { return {a.x*s,a.y*s,a.z*s}; }
inline double dot(Vec a,Vec b) { return a.x*b.x+a.y*b.y+a.z*b.z; }
inline Vec cross(Vec a,Vec b) { return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x}; }
inline double norm(Vec a) { return std::sqrt(dot(a,a)); }
inline Vec unit(Vec a) { return a*(1/norm(a)); }
inline Vec reflected(Vec d,Vec n) { return d-n*(2*dot(d,n)); }
}
