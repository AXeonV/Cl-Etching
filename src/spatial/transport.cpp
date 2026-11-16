#include "transport.hpp"
#include "common/traversal.hpp"
#include "common/periodic.hpp"
#include <algorithm>
namespace etch::spatial {
namespace {
constexpr double infinity=std::numeric_limits<double>::infinity();
}

std::optional<Cell> Transport::reaction_cell(Vec3 point,Vec3 normal,Vec3 origin) {
    int x=static_cast<int>(std::floor(point.x)),y=static_cast<int>(std::floor(point.y)),z=static_cast<int>(std::floor(point.z));
    if(z>=0&&z<grid.geometry.nz&&grid.boundary[grid.index(x,y,z)])return Cell{grid.wrap_x(x),grid.wrap_y(y),z};
    double nearest=infinity;std::optional<Cell> result;
    for(Vec3 d:{normal*-1,normal}) {
        auto interval=ray_box_interval<3>(point,d,origin,surface.fit_radius()+0.5);
        if(!interval)continue;
        GridWalk<3> walk(point,d);
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
        Cell c=grid.coordinates(cell);Vec3 center{c.x+0.5,c.y+0.5,c.z+0.5};
        center.x+=grid.geometry.nx*std::round((point.x-center.x)/grid.geometry.nx);
        center.y+=grid.geometry.ny*std::round((point.y-center.y)/grid.geometry.ny);
        Vec3 offset=center-point;double distance=dot(offset,normal);
        double error=std::max(0.0,dot(offset,offset)-distance*distance);
        if(error<transverse||(error==transverse&&std::abs(distance)<along)) {
            transverse=error;along=std::abs(distance);result=c;
        }
    }
    if(result){++diagnostics.off_normal_mapped_hits;diagnostics.max_mapping_offset=std::max(diagnostics.max_mapping_offset,std::sqrt(transverse));}
    return result;
}
std::optional<Hit> Transport::next_hit(const Particle& particle) {
    Vec3 p=particle.position,d=particle.direction;int n=surface.resolution(),nz=grid.geometry.nz*n;
    if(d.z==0&&!surface.can_hit_horizontal_ray(p.z)){++diagnostics.periodic_orbits;diagnostics.last_end=RayEnd::PeriodicOrbit;return std::nullopt;}
    GridWalk<3> walk(p*n,d*n);double period=d.z==0?horizontal_ray_period(d.x,d.y,grid.geometry.nx,grid.geometry.ny):infinity;
    while(walk.cell[2]>=0&&walk.cell[2]<nz) {
        double start=walk.entry,end=walk.exit();auto cell=walk.cell;
        if(end>start) {
            const auto& field=surface.field(cell[0],cell[1],cell[2]);
            Vec3 origin{double(cell[0]),double(cell[1]),double(cell[2])},a=(p+d*start)*n-origin,b=(p+d*end)*n-origin;
            for(double root:field.intersections({a.x,a.y,a.z},{b.x,b.y,b.z})) {
                double t=start+(end-start)*root;if(t<1e-8)continue;
                Vec3 point=p+d*t,relative=point*n-origin;
                auto gradient=field.gradient({relative.x,relative.y,relative.z});Vec3 n{gradient[0],gradient[1],gradient[2]};
                if(norm(n)==0)continue;
                Vec3 normal=unit(n);if(dot(d,normal)>=-1e-10)continue;auto owner=reaction_cell(point,normal,point);
                if(!owner){++diagnostics.mapping_failures;continue;}
                int mapped=grid.index(owner->x,owner->y,owner->z);
                if(grid.index(int(std::floor(point.x)),int(std::floor(point.y)),int(std::floor(point.z)))!=mapped)++diagnostics.normal_mapped_hits;
                return Hit{point,normal,mapped};
            }
        }
        walk.advance();
        if(walk.entry>=period){++diagnostics.periodic_orbits;diagnostics.last_end=RayEnd::PeriodicOrbit;return std::nullopt;}
    }
    diagnostics.last_end=walk.cell[2]>=nz?RayEnd::MissedSurface:RayEnd::Escaped;
    if(diagnostics.last_end==RayEnd::MissedSurface)++diagnostics.interior_exits;
    return std::nullopt;
}
}
