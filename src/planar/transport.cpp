#include "transport.hpp"
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

std::optional<Cell> Transport::reaction_cell(Vec point, Vec normal, Vec fit_origin) const {
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
    return mapped;
}

std::optional<Hit> Transport::next_hit(const Particle& particle) {
    Vec p=particle.position, d=particle.direction;
    int nx=grid.geometry.nx, nz=grid.geometry.nz;
    GridWalk walk(p,d);
    while (walk.z>=0 && walk.z<nz) {
        double nearest=infinity;
        std::optional<Hit> closest;
        // Enumerate every fitted patch whose support can overlap this ray segment.
        for (int cell : surface.candidates(walk.x,walk.z)) {
            int x=cell%nx, z=cell/nx;
            Curve curve=surface.curve(x,z);
            if (!curve.defined) continue;
            curve.origin.x += nx*std::round((walk.x+0.5-curve.origin.x)/nx);
            auto interval=local_interval(p,d,curve.origin,surface.fit_radius()+0.5);
            if (!interval) continue;
            double start=std::max({1e-8,walk.entry,interval->first});
            double end=std::min({walk.tx,walk.tz,interval->last});
            if (start>end) continue;
            auto distance=curve.intersection(p,d,start,end);
            if (!distance) { ++no_local_intersections; continue; }
            if (*distance>=nearest) continue;
            Vec point=p+d*(*distance);
            Vec normal=curve.normal(dot(point-curve.origin,curve.tangent));
            auto mapped=reaction_cell(point,normal,curve.origin);
            if (!mapped) { ++mapping_failures; continue; }
            nearest=*distance; closest=Hit{point,normal,grid.index(mapped->x,mapped->z)};
        }
        if (closest) {
            if (grid.index(static_cast<int>(std::floor(closest->position.x)),
                           static_cast<int>(std::floor(closest->position.z)))!=closest->cell) ++normal_mapped_hits;
            return closest;
        }
        walk.advance();
        if (d.z==0 && walk.entry>=nx) {
            ++periodic_orbits;last_end=RayEnd::PeriodicOrbit;return std::nullopt;
        }
    }
    last_end=walk.z>=nz?RayEnd::MissedSurface:RayEnd::Escaped;
    if (last_end==RayEnd::MissedSurface) ++interior_exits;
    return std::nullopt;
}
}
