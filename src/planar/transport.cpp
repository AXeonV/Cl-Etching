#include "transport.hpp"
#include "common/implicit.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

namespace etch::planar {
namespace {
constexpr double infinity = std::numeric_limits<double>::infinity();

// Exact DDA cell intervals: advancing grid indices does not skip any ray segment.
struct GridWalk {
    int x, z, sx, sz, corner_step = 0;
    double entry = 0, tx, tz, dx, dz;
    GridWalk(Vec p, Vec d) : x(static_cast<int>(std::floor(p.x))), z(static_cast<int>(std::floor(p.z))),
        sx(d.x > 0 ? 1 : -1), sz(d.z > 0 ? 1 : -1),
        tx(d.x == 0 ? infinity : (x + (sx > 0) - p.x) / d.x),
        tz(d.z == 0 ? infinity : (z + (sz > 0) - p.z) / d.z),
        dx(d.x == 0 ? infinity : 1 / std::abs(d.x)),
        dz(d.z == 0 ? infinity : 1 / std::abs(d.z)) {}
    void advance() {
        // At a grid corner, include both cells touched only at that corner.
        if (corner_step == 1) { x -= sx; z += sz; corner_step = 2; return; }
        if (corner_step == 2) { x += sx; tx += dx; tz += dz; corner_step = 0; return; }
        entry = std::min(tx,tz);
        if (tx == tz) { x += sx; corner_step = 1; }
        else if (tx < tz) { x += sx; tx += dx; }
        else { z += sz; tz += dz; }
    }
};

struct Interval { double first = 0, last = infinity; };
std::optional<Interval> local_interval(Vec p, Vec d, Vec center, double half_width) {
    Interval interval;
    for (auto axis : {std::array<double, 3>{p.x, d.x, center.x}, {p.z, d.z, center.z}}) {
        auto [position, direction, middle] = axis;
        if (direction == 0) {
            if (position < middle - half_width || position > middle + half_width) return std::nullopt;
        } else {
            double first = (middle - half_width - position) / direction;
            double last = (middle + half_width - position) / direction;
            if (first > last) std::swap(first, last);
            interval.first = std::max(interval.first, first);
            interval.last = std::min(interval.last, last);
        }
    }
    if (interval.last < interval.first) return std::nullopt;
    return interval;
}
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
        auto interval = local_interval(point, direction, fit_origin, surface.fit_radius() + 0.5);
        if (!interval) continue;
        GridWalk walk(point, direction);
        while (walk.z >= 0 && walk.z < grid.geometry.nz && walk.entry <= interval->last) {
            if (grid.boundary[grid.index(walk.x, walk.z)] && walk.entry < nearest) {
                nearest = walk.entry;
                mapped = Cell{grid.wrap(walk.x), walk.z};
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
    if(mapped){++off_normal_mapped_hits;max_mapping_offset=std::max(max_mapping_offset,std::sqrt(transverse));}
    return mapped;
}

std::optional<Hit> Transport::next_hit(const Particle& particle) {
    Vec p=particle.position,d=particle.direction;int n=surface.resolution();GridWalk walk(p*n,d*n);int nz=grid.geometry.nz*n;
    while(walk.z>=0&&walk.z<nz) {
        double start=walk.entry,end=std::min(walk.tx,walk.tz);
        if(end>start) {
            const auto& field=surface.field(walk.x,walk.z);
            Vec a=(p+d*start)*n-Vec{double(walk.x),double(walk.z)},b=(p+d*end)*n-Vec{double(walk.x),double(walk.z)};
            for(double root:field.intersections({a.x,a.z},{b.x,b.z})) {
                double t=start+(end-start)*root;if(t<1e-8)continue;
                Vec point=p+d*t,relative=point*n-Vec{double(walk.x),double(walk.z)};
                auto gradient=field.gradient({relative.x,relative.z});Vec n{gradient[0],gradient[1]};
                if(norm(n)==0)continue;
                Vec normal=unit(n);if(dot(d,normal)>=-1e-10)continue;auto owner=reaction_cell(point,normal,point);
                if(!owner){++mapping_failures;continue;}
                int cell=grid.index(owner->x,owner->z);
                if(grid.index(int(std::floor(point.x)),int(std::floor(point.z)))!=cell)++normal_mapped_hits;
                return Hit{point,normal,cell};
            }
        }
        walk.advance();
        if(d.z==0&&walk.entry>=grid.geometry.nx){++periodic_orbits;last_end=RayEnd::PeriodicOrbit;return std::nullopt;}
    }
    last_end=walk.z>=nz?RayEnd::MissedSurface:RayEnd::Escaped;
    if(last_end==RayEnd::MissedSurface)++interior_exits;
    return std::nullopt;
}
}
