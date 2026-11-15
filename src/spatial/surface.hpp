#pragma once
#include "grid.hpp"
#include "common/implicit.hpp"
#include <unordered_map>
namespace etch::spatial {
struct Projection { double u,v,distance; };
struct Patch {
    Vec origin,u_axis,v_axis,outward;
    // w = a*u^2 + b*v^2 + c*u*v + d*u + e*v + f.
    double a=0,b=0,c=0,d=0,e=0,f=0;
    bool defined=true;
    double height(double u,double v) const { return a*u*u+b*v*v+c*u*v+d*u+e*v+f; }
    Vec position(double u,double v) const { return origin+u_axis*u+v_axis*v+outward*height(u,v); }
    Vec normal(double u,double v) const;
    Projection project(Vec point) const;
};
class Surface {
    Grid& grid;
    int radius,subdivisions;
    using Field=ImplicitCell<3>;
    std::vector<Field::Node> node_values;
    std::unordered_map<int,Field> fields;
    std::vector<uint8_t> dirty_nodes;
    std::unordered_map<int,Patch> patches;
    std::vector<std::vector<int>> nearby_cells;
    std::vector<uint8_t> dirty_candidates;
public:
    long long fit_count=0,rank_deficient_fits=0,undefined_normal_fits=0;
    Surface(Grid& grid,int radius=3,int subdivisions=1);
    Field::Node node(int x,int y,int z);
    const Field& field(int x,int y,int z);
    int resolution() const { return subdivisions; }
    int fit_radius() const { return radius; }
    const Patch& patch(int x,int y,int z);
    const std::vector<int>& candidates(int x,int y,int z);
    void invalidate(int x,int y,int z);
    bool can_hit_horizontal_ray(double z);
};
}
