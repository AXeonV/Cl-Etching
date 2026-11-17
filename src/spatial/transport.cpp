#include "transport.hpp"
#include "common/traversal.hpp"
#include "common/periodic.hpp"
#include <algorithm>

namespace etch::spatial {
namespace {
constexpr double infinity = std::numeric_limits<double>::infinity();
}

std::optional<Cell> Transport::reaction_cell(Vec3 point, Vec3 normal, Vec3 origin) {
  int x = static_cast<int>(std::floor(point.x)), y = static_cast<int>(std::floor(point.y)),
      z = static_cast<int>(std::floor(point.z));
  if (z >= 0 && z < grid.geometry.nz && grid.boundary[grid.index(x, y, z)])
    return Cell{grid.wrap_x(x), grid.wrap_y(y), z};

  // 优先沿碰撞法线两侧寻找拟合邻域内最近的表面格。
  double nearest = infinity;
  std::optional<Cell> result;
  for (Vec3 d : {normal * -1, normal}) {
    auto interval = ray_box_interval<3>(point, d, origin, surface.fit_radius() + 0.5);
    if (!interval)
      continue;
    GridWalk<3> walk(point, d);
    while (walk.cell[2] >= 0 && walk.cell[2] < grid.geometry.nz && walk.entry <= interval->last) {
      auto p = walk.cell;
      if (grid.boundary[grid.index(p[0], p[1], p[2])] && walk.entry < nearest) {
        nearest = walk.entry;
        result = Cell{grid.wrap_x(p[0]), grid.wrap_y(p[1]), p[2]};
        break;
      }
      walk.advance();
    }
  }
  if (result)
    return result;

  // 法线未穿过边界格时，按垂距、沿法线距离依次完成离散映射。
  double transverse = infinity, along = infinity;
  for (int cell : surface.candidates(x, y, z)) {
    Cell c = grid.coordinates(cell);
    Vec3 center{c.x + 0.5, c.y + 0.5, c.z + 0.5};
    center.x += grid.geometry.nx * std::round((point.x - center.x) / grid.geometry.nx);
    center.y += grid.geometry.ny * std::round((point.y - center.y) / grid.geometry.ny);
    Vec3 offset = center - point;
    double distance = dot(offset, normal);
    double error = std::max(0.0, dot(offset, offset) - distance * distance);
    if (error < transverse || (error == transverse && std::abs(distance) < along)) {
      transverse = error;
      along = std::abs(distance);
      result = c;
    }
  }
  if (result) {
    ++diagnostics.off_normal_mapped_hits;
    diagnostics.max_mapping_offset =
      std::max(diagnostics.max_mapping_offset, std::sqrt(transverse));
  }
  return result;
}

// 在连续碰撞场中求交，再映射到离散反应格；两步使用同一个交点。
std::optional<Hit> Transport::next_hit(const Particle& particle) {
  Vec3 ray_origin = particle.position, direction = particle.direction;
  int resolution = surface.resolution(), reconstruction_layers = grid.geometry.nz * resolution;
  if (direction.z == 0 && !surface.can_hit_horizontal_ray(ray_origin.z)) {
    ++diagnostics.periodic_orbits;
    diagnostics.last_end = RayEnd::PeriodicOrbit;
    return std::nullopt;
  }
  GridWalk<3> walk(ray_origin * resolution, direction * resolution);
  double period = direction.z == 0 ? horizontal_ray_period(
                                       direction.x, direction.y, grid.geometry.nx, grid.geometry.ny)
                                   : infinity;
  while (walk.cell[2] >= 0 && walk.cell[2] < reconstruction_layers) {
    double interval_start = walk.entry, interval_end = walk.exit();
    auto reconstruction_cell = walk.cell;
    if (interval_end > interval_start) {
      const auto& field =
        surface.field(reconstruction_cell[0], reconstruction_cell[1], reconstruction_cell[2]);
      Vec3 cell_origin{double(reconstruction_cell[0]), double(reconstruction_cell[1]),
        double(reconstruction_cell[2])},
        segment_start = (ray_origin + direction * interval_start) * resolution - cell_origin,
        segment_end = (ray_origin + direction * interval_end) * resolution - cell_origin;

      // 每个 DDA 区间保留全部交点，仅接受沿运动方向进入固体的交点。
      for (double root : field.intersections({segment_start.x, segment_start.y, segment_start.z},
             {segment_end.x, segment_end.y, segment_end.z})) {
        double travel = interval_start + (interval_end - interval_start) * root;
        // 反射从原交点继续，排除零程自撞。
        if (travel < 1e-8)
          continue;
        Vec3 point = ray_origin + direction * travel,
             local_point = point * resolution - cell_origin;
        auto gradient_components = field.gradient({local_point.x, local_point.y, local_point.z});
        Vec3 field_gradient{gradient_components[0], gradient_components[1], gradient_components[2]};
        if (norm(field_gradient) == 0)
          continue;
        Vec3 normal = unit(field_gradient);
        if (dot(direction, normal) >= -1e-10)
          continue;

        auto surface_cell = reaction_cell(point, normal, point);
        if (!surface_cell) {
          ++diagnostics.mapping_failures;
          continue;
        }
        int material_cell = grid.index(surface_cell->x, surface_cell->y, surface_cell->z);
        if (grid.index(int(std::floor(point.x)), int(std::floor(point.y)),
              int(std::floor(point.z))) != material_cell)
          ++diagnostics.normal_mapped_hits;
        return Hit{point, normal, material_cell};
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
    walk.cell[2] >= reconstruction_layers ? RayEnd::MissedSurface : RayEnd::Escaped;
  if (diagnostics.last_end == RayEnd::MissedSurface)
    ++diagnostics.interior_exits;
  return std::nullopt;
}
}
