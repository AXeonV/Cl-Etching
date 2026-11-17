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

  // 交点不属于表面格时，沿法线两侧搜索拟合邻域内最近的表面格（刘婉正文第 38 页）。
  double nearest = infinity;
  std::optional<Cell> mapped;
  for (Vec direction : {normal * -1, normal}) {
    auto interval = ray_box_interval<2>(point, direction, fit_origin, surface.fit_radius() + 0.5);
    if (!interval)
      continue;
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
  if (mapped)
    return mapped;

  // 法线可能不穿过任何离散表面格；此时按格中心到法线的垂距选择反应格。
  double transverse = infinity, along = infinity;
  for (int cell : surface.candidates(x, z)) {
    Vec center{cell % grid.geometry.nx + 0.5, cell / grid.geometry.nx + 0.5};
    center.x += grid.geometry.nx * std::round((point.x - center.x) / grid.geometry.nx);
    Vec offset = center - point;
    double distance = dot(offset, normal);
    double error = std::max(0.0, dot(offset, offset) - distance * distance);
    if (error < transverse || (error == transverse && std::abs(distance) < along)) {
      transverse = error;
      along = std::abs(distance);
      mapped = Cell{cell % grid.geometry.nx, cell / grid.geometry.nx};
    }
  }
  if (mapped) {
    ++diagnostics.off_normal_mapped_hits;
    diagnostics.max_mapping_offset =
      std::max(diagnostics.max_mapping_offset, std::sqrt(transverse));
  }
  return mapped;
}

// 在连续碰撞场中求交，再映射到离散反应格；两步使用同一个交点。
std::optional<Hit> Transport::next_hit(const Particle& particle) {
  Vec ray_origin = particle.position, direction = projection(particle.direction);
  if (direction.x == 0 && direction.z == 0) {
    ++diagnostics.periodic_orbits;
    diagnostics.last_end = RayEnd::PeriodicOrbit;
    return std::nullopt;
  }
  double period = direction.z == 0 ? grid.geometry.nx / std::abs(direction.x) : infinity;
  double projected_speed = norm(direction);
  int resolution = surface.resolution();
  GridWalk<2> walk(ray_origin * resolution, direction * resolution);
  int reconstruction_layers = grid.geometry.nz * resolution;
  while (walk.cell[1] >= 0 && walk.cell[1] < reconstruction_layers) {
    double interval_start = walk.entry, interval_end = walk.exit();
    if (interval_end > interval_start) {
      const auto& field = surface.field(walk.cell[0], walk.cell[1]);
      Vec segment_start = (ray_origin + direction * interval_start) * resolution -
                          Vec{double(walk.cell[0]), double(walk.cell[1])},
          segment_end = (ray_origin + direction * interval_end) * resolution -
                        Vec{double(walk.cell[0]), double(walk.cell[1])};

      // 每个 DDA 区间保留全部交点，仅接受沿运动方向进入固体的交点。
      for (double root :
        field.intersections({segment_start.x, segment_start.z}, {segment_end.x, segment_end.z})) {
        double travel = interval_start + (interval_end - interval_start) * root;
        // 反射从原交点继续，排除零程自撞。
        if (travel < 1e-8)
          continue;
        Vec point = ray_origin + direction * travel,
            local_point = point * resolution - Vec{double(walk.cell[0]), double(walk.cell[1])};
        auto gradient_components = field.gradient({local_point.x, local_point.z});
        Vec field_gradient{gradient_components[0], gradient_components[1]};
        if (norm(field_gradient) == 0)
          continue;
        Vec normal = unit(field_gradient);
        if (dot(direction, normal) >= -1e-10 * projected_speed)
          continue;

        auto surface_cell = reaction_cell(point, normal, point);
        if (!surface_cell) {
          ++diagnostics.mapping_failures;
          continue;
        }
        int material_cell = grid.index(surface_cell->x, surface_cell->z);
        if (grid.index(int(std::floor(point.x)), int(std::floor(point.z))) != material_cell)
          ++diagnostics.normal_mapped_hits;
        return Hit{point, {normal.x, 0, normal.z}, material_cell};
      }
    }

    walk.advance();

    // 走完无碰撞的横向周期后不会遇到新表面，结束该周期轨迹。
    if (walk.entry >= period) {
      ++diagnostics.periodic_orbits;
      diagnostics.last_end = RayEnd::PeriodicOrbit;
      return std::nullopt;
    }
  }
  diagnostics.last_end =
    walk.cell[1] >= reconstruction_layers ? RayEnd::MissedSurface : RayEnd::Escaped;
  if (diagnostics.last_end == RayEnd::MissedSurface)
    ++diagnostics.interior_exits;
  return std::nullopt;
}
}
