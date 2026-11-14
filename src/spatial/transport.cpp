#include "transport.hpp"
#include "traversal.hpp"
#include <algorithm>
namespace etch::spatial {
namespace {
constexpr double infinity=std::numeric_limits<double>::infinity();
struct Interval { double first=0,last=infinity; };
std::optional<Interval> local_interval(Vec p,Vec d,Vec center,double half) {
    Interval result;
    for(int axis=0;axis<3;++axis) {
        if(d[axis]==0) { if(p[axis]<center[axis]-half||p[axis]>center[axis]+half)return std::nullopt; }
        else {
            double first=(center[axis]-half-p[axis])/d[axis],last=(center[axis]+half-p[axis])/d[axis];
            if(first>last)std::swap(first,last);
            result.first=std::max(result.first,first);result.last=std::min(result.last,last);
        }
    }
    if(result.last<result.first)return std::nullopt;
    return result;
}
double horizontal_period(Vec d,const Geometry& g) {
    if(d.z!=0)return infinity;
    if(d.x==0)return g.ny/std::abs(d.y);
    if(d.y==0)return g.nx/std::abs(d.x);
    double x=g.nx/std::abs(d.x),y=g.ny/std::abs(d.y);
    return x==y?x:infinity;
}
}
std::optional<Cell> Transport::reaction_cell(Vec point,Vec normal,Vec origin) const {
    int x=static_cast<int>(std::floor(point.x)),y=static_cast<int>(std::floor(point.y)),z=static_cast<int>(std::floor(point.z));
    if(z>=0&&z<grid.geometry.nz&&grid.boundary[grid.index(x,y,z)])return Cell{grid.wrap_x(x),grid.wrap_y(y),z};
    double nearest=infinity;std::optional<Cell> result;
    for(Vec d:{normal*-1,normal}) {
        auto interval=local_interval(point,d,origin,surface.fit_radius()+0.5);
        if(!interval)continue;
        VoxelWalk walk(point,d);
        while(walk.cell[2]>=0&&walk.cell[2]<grid.geometry.nz&&walk.entry<=interval->last) {
            auto p=walk.cell;
            if(grid.boundary[grid.index(p[0],p[1],p[2])]&&walk.entry<nearest) {
                nearest=walk.entry;result=Cell{grid.wrap_x(p[0]),grid.wrap_y(p[1]),p[2]};break;
            }
            walk.advance();
        }
    }
    return result;
}
std::optional<Hit> Transport::next_hit(const Particle& particle) {
    Vec p=particle.position,d=particle.direction;
    int nx=grid.geometry.nx,ny=grid.geometry.ny,nz=grid.geometry.nz;
    if(d.z==0 && !surface.can_hit_horizontal_ray(p.z,d)) {
        ++periodic_orbits;last_end=RayEnd::PeriodicOrbit;return std::nullopt;
    }
    VoxelWalk walk(p,d);double period=horizontal_period(d,grid.geometry);
    while(walk.cell[2]>=0&&walk.cell[2]<nz) {
        double nearest=infinity;std::optional<Hit> best;
        auto coordinate=walk.cell;
        for(int cell:surface.candidates(coordinate[0],coordinate[1],coordinate[2])) {
            Cell index=grid.coordinates(cell);Patch patch=surface.patch(index.x,index.y,index.z);
            if(!patch.defined)continue;
            patch.origin.x+=nx*std::round((coordinate[0]+0.5-patch.origin.x)/nx);
            patch.origin.y+=ny*std::round((coordinate[1]+0.5-patch.origin.y)/ny);
            auto interval=local_interval(p,d,patch.origin,surface.fit_radius()+0.5);
            if(!interval)continue;
            double start=std::max({1e-8,walk.entry,interval->first}),end=std::min(walk.exit(),interval->last);
            if(start>end)continue;
            auto distance=patch.intersection(p,d,start,end);
            if(!distance){++no_local_intersections;continue;}
            if(*distance>=nearest)continue;
            Vec point=p+d*(*distance),relative=point-patch.origin;
            Vec normal=patch.normal(dot(relative,patch.u_axis),dot(relative,patch.v_axis));
            auto mapped=reaction_cell(point,normal,patch.origin);
            if(!mapped){++mapping_failures;continue;}
            nearest=*distance;best=Hit{point,normal,grid.index(mapped->x,mapped->y,mapped->z)};
        }
        if(best) {
            if(grid.index(static_cast<int>(std::floor(best->position.x)),static_cast<int>(std::floor(best->position.y)),
                          static_cast<int>(std::floor(best->position.z)))!=best->cell)++normal_mapped_hits;
            return best;
        }
        walk.advance();
        if(walk.entry>=period){++periodic_orbits;last_end=RayEnd::PeriodicOrbit;return std::nullopt;}
    }
    last_end=walk.cell[2]>=nz?RayEnd::MissedSurface:RayEnd::Escaped;
    if(last_end==RayEnd::MissedSurface)++interior_exits;
    return std::nullopt;
}
}
