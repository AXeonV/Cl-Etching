#include "transport.hpp"
#include "common/traversal.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

namespace etch::planar {
namespace {
constexpr double infinity = std::numeric_limits<double>::infinity();

}

std::optional<Cell> Transport::reaction_cell(Vec point, Vec normal, Vec fit_origin) {
    int x = static_cast<int>(std::floor(point.x)), z = static_cast<int>(std::floor(point.z));
    if (z >= 0 && z < grid.geometry.nz && grid.boundary[grid.index(x, z)])
        return Cell{grid.wrap(x), z};

    // If the intersection is not in a surface cell, use the closest surface cell
    // intersected by its normal line within the fitted neighborhood (Liu, p.38).
    double nearest = infinity;
    std::optional<Cell> mapped;
    for (Vec direction : {normal * -1, normal}) {
        auto interval = ray_box_interval<2>(point, direction, fit_origin, surface.fit_radius() + 0.5);
        if (!interval) continue;
        GridWalk<2> walk(point, direction);
        while (walk.cell[1] >= 0 && walk.cell[1] < grid.geometry.nz && walk.entry <= interval->last) {
            if (grid.boundary[grid.index(walk.cell[0], walk.cell[1])] && walk.entry < nearest) {
                nearest = walk.entry;
                mapped = Cell{grid.wrap(walk.cell[0]), walk.cell[1]};
                break;
            }
            walk.advance();
        }
    }
    if(mapped)return mapped;
    // A fitted normal need not pass through a boundary voxel. Complete the discrete
    // projection by choosing the closest surface-cell center to the normal line.
    double transverse=infinity,along=infinity;
    for(int cell:surface.candidates(x,z)) {
        Vec center{cell%grid.geometry.nx+0.5,cell/grid.geometry.nx+0.5};
        center.x+=grid.geometry.nx*std::round((point.x-center.x)/grid.geometry.nx);
        Vec offset=center-point;double distance=dot(offset,normal);
        double error=std::max(0.0,dot(offset,offset)-distance*distance);
        if(error<transverse||(error==transverse&&std::abs(distance)<along)) {
            transverse=error;along=std::abs(distance);mapped=Cell{cell%grid.geometry.nx,cell/grid.geometry.nx};
        }
    }
    if(mapped){++diagnostics.off_normal_mapped_hits;diagnostics.max_mapping_offset=std::max(diagnostics.max_mapping_offset,std::sqrt(transverse));}
    return mapped;
}

std::optional<Hit> Transport::next_hit(const Particle& particle) {
    Vec p=particle.position,d=projection(particle.direction);
    if(d.x==0&&d.z==0){++diagnostics.periodic_orbits;diagnostics.last_end=RayEnd::PeriodicOrbit;return std::nullopt;}
    double period=d.z==0?grid.geometry.nx/std::abs(d.x):infinity;
    double projected_speed=norm(d);
    int n=surface.resolution();GridWalk<2> walk(p*n,d*n);int nz=grid.geometry.nz*n;
    while(walk.cell[1]>=0&&walk.cell[1]<nz) {
        double start=walk.entry,end=walk.exit();
        if(end>start) {
            const auto& field=surface.field(walk.cell[0],walk.cell[1]);
            Vec a=(p+d*start)*n-Vec{double(walk.cell[0]),double(walk.cell[1])},b=(p+d*end)*n-Vec{double(walk.cell[0]),double(walk.cell[1])};
            for(double root:field.intersections({a.x,a.z},{b.x,b.z})) {
                double t=start+(end-start)*root;if(t<1e-8)continue;
                Vec point=p+d*t,relative=point*n-Vec{double(walk.cell[0]),double(walk.cell[1])};
                auto gradient=field.gradient({relative.x,relative.z});Vec n{gradient[0],gradient[1]};
                if(norm(n)==0)continue;
                Vec normal=unit(n);if(dot(d,normal)>=-1e-10*projected_speed)continue;auto owner=reaction_cell(point,normal,point);
                if(!owner){++diagnostics.mapping_failures;continue;}
                int cell=grid.index(owner->x,owner->z);
                if(grid.index(int(std::floor(point.x)),int(std::floor(point.z)))!=cell)++diagnostics.normal_mapped_hits;
                return Hit{point,{normal.x,0,normal.z},cell};
            }
        }
        walk.advance();
        if(walk.entry>=period){++diagnostics.periodic_orbits;diagnostics.last_end=RayEnd::PeriodicOrbit;return std::nullopt;}
    }
    diagnostics.last_end=walk.cell[1]>=nz?RayEnd::MissedSurface:RayEnd::Escaped;
    if(diagnostics.last_end==RayEnd::MissedSurface)++diagnostics.interior_exits;
    return std::nullopt;
}
}
