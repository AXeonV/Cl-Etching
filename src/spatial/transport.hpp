#pragma once
#include <optional>
#include "source.hpp"
#include "surface.hpp"
#include "common/events.hpp"

namespace etch::spatial {
using Hit = SurfaceHit<Vec3>;

class Transport {
  Grid& grid;
  Surface& surface;

public:
  TrackingDiagnostics diagnostics;

  Transport(Grid& grid, Surface& surface) : grid(grid), surface(surface) {}

  std::optional<Hit> next_hit(const Particle& particle);
  std::optional<Cell> reaction_cell(Vec3 point, Vec3 normal, Vec3 fit_origin);
};
}
