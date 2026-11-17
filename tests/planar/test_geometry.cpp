#include "planar/backend.hpp"
#include <iostream>
#include <algorithm>
#include <limits>

using namespace etch;
using namespace etch::planar;

// 独立参考解：对共享场采样，二分求解首次由气体进入固体的变号区间。
std::optional<Hit> sampled_hit(
  Grid& grid, Surface& surface, Transport& transport, Particle ray, double limit) {
  auto evaluate = [&](double distance) {
    Vec p = (ray.position + projection(ray.direction) * distance) * surface.resolution();
    int x = int(std::floor(p.x)), z = int(std::floor(p.z));
    return surface.field(x, z).value({p.x - x, p.z - z});
  };
  double first = 1e-7, previous = evaluate(first);
  for (double last = std::min(first + 1.0 / 64, limit); first < limit;
    last = std::min(last + 1.0 / 64, limit)) {
    double current = evaluate(last);
    if (previous > 0 && current <= 0) {
      double left = first, right = last;
      for (int iteration = 0; iteration < 50; ++iteration) {
        double middle = (left + right) / 2;
        if (evaluate(middle) > 0)
          left = middle;
        else
          right = middle;
      }
      Vec point = ray.position + projection(ray.direction) * ((left + right) / 2),
          p = point * surface.resolution();
      int x = int(std::floor(p.x)), z = int(std::floor(p.z));
      auto g = surface.field(x, z).gradient({p.x - x, p.z - z});
      Vec normal = unit(Vec{g[0], g[1]});
      auto owner = transport.reaction_cell(point, normal, point);
      if (owner)
        return Hit{point, {normal.x, 0, normal.z}, grid.index(owner->x, owner->z)};
      return std::nullopt;
    }
    first = last;
    previous = current;
  }
  return std::nullopt;
}

int mask_corner_reflections() {
  Grid grid(Geometry{256, 512, 60, 24, 72});
  Surface surface(grid);
  Transport transport(grid, surface);
  int misses = 0, self_hits = 0, reference_errors = 0;
  for (int i = 0; i < 80; ++i)
    for (int j = 0; j < 80; ++j) {
      double z = 55 + (i + 0.5) * 6 / 80, dx = 0.01 + (j + 0.5) * 0.3 / 80;
      Vec3 direction = unit(Vec3{dx, 0, 1});
      Particle ray{{164.5 - (z - 0.25) * direction.x / direction.z, 0.25}, direction, 100};
      auto first = transport.next_hit(ray);
      if (!first) {
        ++misses;
        continue;
      }
      ray.position = first->position;
      ray.direction = reflected(ray.direction, first->normal);
      auto second = transport.next_hit(ray);
      if (!second) {
        if (transport.diagnostics.last_end == RayEnd::MissedSurface)
          ++misses;
      } else if (norm(second->position - first->position) < 1e-4)
        ++self_hits;
      if (i % 8 == 0 && j % 8 == 0) {
        auto reference = sampled_hit(
          grid, surface, transport, ray, second ? norm(second->position - ray.position) + 0.01 : 0);
        if (bool(reference) != bool(second) ||
            (reference && norm(reference->position - second->position) > 1e-6))
          ++reference_errors;
      }
    }
  std::cout << (misses || self_hits || reference_errors ? "FAIL " : "PASS ")
            << "6400 mask-corner reflection rays: misses=" << misses << ", self_hits=" << self_hits
            << ", sampled-field-reference errors=" << reference_errors << "/100\n";
  return misses || self_hits || reference_errors;
}

int main() {
  // 几何诊断算例，输入独立于论文图中的工况。
  Geometry geometry{256, 512, 60, 24, 72};
  Grid flat(geometry);
  Surface surface(flat);
  Transport transport(flat, surface);
  int wrong_cells = 0;
  for (int i = 0; i < 100; ++i) {
    double x = 125 + (i + 0.5) / 100;
    auto hit = transport.next_hit({{x, 1}, {0, 0, 1}, 100});
    if (!hit || hit->cell % geometry.nx != static_cast<int>(std::floor(x)))
      ++wrong_cells;
  }
  std::cout << (wrong_cells ? "FAIL " : "PASS ") << "flat-surface hit ownership: " << wrong_cells
            << "/100 assigned to a different cell\n";

  Settings settings{geometry};
  settings.max_bounces = 500;
  double sigma = 10 * std::acos(-1.0) / 180 / 2.354820045;
  Source source{gaussian_direction(sigma), constant(100)};
  Simulation simulation(settings, source, std::make_unique<YieldReaction>());
  simulation.advance(100000);
  std::cout << (simulation.transport.diagnostics.interior_exits ? "FAIL " : "PASS ")
            << "surface interception: " << simulation.transport.diagnostics.interior_exits
            << "/100000 rays exited through the substrate\n";
  int corner_errors = mask_corner_reflections();
  return wrong_cells || simulation.transport.diagnostics.interior_exits || corner_errors ? 1 : 0;
}
