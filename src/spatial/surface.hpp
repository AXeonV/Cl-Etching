#pragma once
#include "grid.hpp"
#include <optional>
#include <unordered_map>
#include <limits>
namespace etch::spatial {
struct Projection { double u,v,distance; };
struct Patch {
    Vec origin,u_axis,v_axis,outward;
    // w = a*u^2 + b*v^2 + c*u*v + d*u + e*v + f.
    double a=0,b=0,c=0,d=0,e=0,f=0;
    double u_min=-std::numeric_limits<double>::infinity(),u_max=std::numeric_limits<double>::infinity();
    double v_min=-std::numeric_limits<double>::infinity(),v_max=std::numeric_limits<double>::infinity();
    bool defined=true;
    double height(double u,double v) const { return a*u*u+b*v*v+c*u*v+d*u+e*v+f; }
    Vec position(double u,double v) const { return origin+u_axis*u+v_axis*v+outward*height(u,v); }
    Vec normal(double u,double v) const;
    Projection project(Vec point) const;
    std::pair<double,double> z_bounds() const;
    std::optional<double> intersection(Vec p,Vec direction,double start,double end) const;
};
class Surface {
    Grid& grid;
    int radius;
    std::unordered_map<int,Patch> patches;
    std::vector<std::vector<int>> nearby_cells;
    std::vector<uint8_t> dirty_candidates;
public:
    long long fit_count=0,rank_deficient_fits=0,undefined_normal_fits=0;
    Surface(Grid& grid,int radius=3);
    int fit_radius() const { return radius; }
    const Patch& patch(int x,int y,int z);
    const std::vector<int>& candidates(int x,int y,int z);
    void invalidate(int x,int y,int z);
    bool can_hit_horizontal_ray(double z,Vec direction);
};
}
