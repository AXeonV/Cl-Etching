#include "transport.hpp"
#include "common/implicit.hpp"
#include "traversal.hpp"
#include "common/periodic.hpp"
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

}
std::optional<Cell> Transport::reaction_cell(Vec point,Vec normal,Vec origin) {
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
    if(result)return result;
    double transverse=infinity,along=infinity;
    for(int cell:surface.candidates(x,y,z)) {
        Cell c=grid.coordinates(cell);Vec center{c.x+0.5,c.y+0.5,c.z+0.5};
        center.x+=grid.geometry.nx*std::round((point.x-center.x)/grid.geometry.nx);
        center.y+=grid.geometry.ny*std::round((point.y-center.y)/grid.geometry.ny);
        Vec offset=center-point;double distance=dot(offset,normal);
        double error=std::max(0.0,dot(offset,offset)-distance*distance);
        if(error<transverse||(error==transverse&&std::abs(distance)<along)) {
            transverse=error;along=std::abs(distance);result=c;
        }
    }
    if(result){++off_normal_mapped_hits;max_mapping_offset=std::max(max_mapping_offset,std::sqrt(transverse));}
    return result;
}
std::optional<Hit> Transport::next_hit(const Particle& particle) {
    Vec p=particle.position,d=particle.direction;int n=surface.resolution(),nz=grid.geometry.nz*n;
    if(d.z==0&&!surface.can_hit_horizontal_ray(p.z)){++periodic_orbits;last_end=RayEnd::PeriodicOrbit;return std::nullopt;}
    VoxelWalk walk(p*n,d*n);double period=d.z==0?horizontal_ray_period(d.x,d.y,grid.geometry.nx,grid.geometry.ny):infinity;
    while(walk.cell[2]>=0&&walk.cell[2]<nz) {
        double start=walk.entry,end=walk.exit();auto cell=walk.cell;
        if(end>start) {
            const auto& field=surface.field(cell[0],cell[1],cell[2]);
            Vec origin{double(cell[0]),double(cell[1]),double(cell[2])},a=(p+d*start)*n-origin,b=(p+d*end)*n-origin;
            for(double root:field.intersections({a.x,a.y,a.z},{b.x,b.y,b.z})) {
                double t=start+(end-start)*root;if(t<1e-8)continue;
                Vec point=p+d*t,relative=point*n-origin;
                auto gradient=field.gradient({relative.x,relative.y,relative.z});Vec n{gradient[0],gradient[1],gradient[2]};
                if(norm(n)==0)continue;
                Vec normal=unit(n);if(dot(d,normal)>=-1e-10)continue;auto owner=reaction_cell(point,normal,point);
                if(!owner){++mapping_failures;continue;}
                int mapped=grid.index(owner->x,owner->y,owner->z);
                if(grid.index(int(std::floor(point.x)),int(std::floor(point.y)),int(std::floor(point.z)))!=mapped)++normal_mapped_hits;
                return Hit{point,normal,mapped};
            }
        }
        walk.advance();
        if(walk.entry>=period){++periodic_orbits;last_end=RayEnd::PeriodicOrbit;return std::nullopt;}
    }
    last_end=walk.cell[2]>=nz?RayEnd::MissedSurface:RayEnd::Escaped;
    if(last_end==RayEnd::MissedSurface)++interior_exits;
    return std::nullopt;
}
}
