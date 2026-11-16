#pragma once
#include <optional>
#include "source.hpp"
#include "common/events.hpp"
#include "surface.hpp"

namespace etch::planar {
using Hit=SurfaceHit<Vec>;
struct Cell { int x, z; };

class Transport {
    Grid& grid;
    Surface& surface;
public:
    TrackingDiagnostics diagnostics;
    Transport(Grid& grid, Surface& surface) : grid(grid), surface(surface) {}
    std::optional<Hit> next_hit(const Particle& particle);
    std::optional<Cell> reaction_cell(Vec point, Vec normal, Vec fit_origin);
};
}

