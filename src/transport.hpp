#pragma once
#include "source.hpp"
#include "surface.hpp"

namespace etch {
struct Hit { Vec position, normal; int x, z; };
struct Cell { int x, z; };
enum class RayEnd { Escaped, MissedSurface, PeriodicOrbit };
class Transport {
    Grid& grid;
    Surface& surface;
public:
    long long interior_exits = 0;
    long long periodic_orbits = 0;
    long long normal_mapped_hits = 0, mapping_failures = 0, no_local_intersections = 0;
    RayEnd last_end = RayEnd::Escaped;
    Transport(Grid& grid, Surface& surface) : grid(grid), surface(surface) {}
    std::optional<Hit> next_hit(const Particle& particle);
    std::optional<Cell> reaction_cell(Vec point, Vec normal, Vec fit_origin) const;
};
}

