#pragma once
#include "grid.hpp"
#include "source.hpp"
#include "surface.hpp"
#include "transport.hpp"
#include "cell_update.hpp"
#include "common/engine.hpp"

namespace etch::planar {
struct Backend {
  using Geometry = planar::Geometry;
  using Grid = planar::Grid;
  using Surface = planar::Surface;
  using Transport = planar::Transport;
  using CellUpdater = planar::CellUpdater;
  using Source = planar::Source;
  static void write_surface(Surface& surface, const Grid& grid, const std::string& path);
  static void write_field(Surface& surface, const Grid& grid, const std::string& path);
  static void write_geometry(Summary& out, const Geometry& geometry);
};

using Settings = RunSettings<Geometry>;
using Simulation = Engine<Backend>;
}
