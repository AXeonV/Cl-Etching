#pragma once
#include <cmath>
namespace etch {
struct Vec3 {
    double x=0,y=0,z=0;
    double operator[](int axis) const { return axis==0?x:axis==1?y:z; }
};
inline Vec3 operator+(Vec3 a,Vec3 b) { return {a.x+b.x,a.y+b.y,a.z+b.z}; }
inline Vec3 operator-(Vec3 a,Vec3 b) { return {a.x-b.x,a.y-b.y,a.z-b.z}; }
inline Vec3 operator*(Vec3 a,double s) { return {a.x*s,a.y*s,a.z*s}; }
inline double dot(Vec3 a,Vec3 b) { return a.x*b.x+a.y*b.y+a.z*b.z; }
inline Vec3 cross(Vec3 a,Vec3 b) { return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x}; }
inline double norm(Vec3 a) { return std::sqrt(dot(a,a)); }
inline Vec3 unit(Vec3 a) { return a*(1/norm(a)); }
inline Vec3 reflected(Vec3 d,Vec3 n) { return d-n*(2*dot(d,n)); }
}
