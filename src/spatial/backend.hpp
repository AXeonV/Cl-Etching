#pragma once
#include "grid.hpp"
#include "source.hpp"
#include "surface.hpp"
#include "transport.hpp"
#include "cell_update.hpp"
#include "common/engine.hpp"

namespace etch::spatial {
struct Backend {
  using Geometry = spatial::Geometry;
  using Grid = spatial::Grid;
  using Surface = spatial::Surface;
  using Transport = spatial::Transport;
  using CellUpdater = spatial::CellUpdater;
  using Source = spatial::Source;
  static void write_surface(Surface& surface, const Grid& grid, const std::string& path);
  static void write_field(Surface& surface, const Grid& grid, const std::string& path);
  static void write_geometry(Summary& out, const Geometry& geometry);
};

using Settings = RunSettings<Geometry>;
using Simulation = Engine<Backend>;
}
