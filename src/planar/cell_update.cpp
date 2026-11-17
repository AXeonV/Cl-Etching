#include "cell_update.hpp"
#include <algorithm>

namespace etch::planar {
CellUpdater::CellUpdater(Grid& g, Surface& s, double weight, unsigned long long seed)
    : grid(g), surface(s), normal_weight(weight), update(g, diagnostics, seed) {}

void CellUpdater::erase(int cell) {
  int x = cell % grid.geometry.nx, z = cell / grid.geometry.nx;
  grid.remove_empty_cell(x, z);
  surface.invalidate(x, z);
}

CellRemoval CellUpdater::remove(int original, Vec3 collision_normal, Vec collision_position) {
  auto erase_cell = [&](int cell) { erase(cell); };
  // 对称拟合无定义不等于孤立；只提前释放真正没有面连接的格子。
  if (update.isolated(original))
    return update.detach(original, erase_cell);
  int x = original % grid.geometry.nx, z = original / grid.geometry.nx;
  Curve curve = surface.curve(x, z);
  Vec collision{collision_normal.x, collision_normal.z}, normal = collision;
  bool use_curve = curve.defined;
  double orientation = 1;
  if (!curve.defined)
    ++diagnostics.undefined_fit_removals;

  // 先取碰撞点附近的周期像，再按实际接触位置校准拟合的外向符号。
  if (use_curve) {
    collision_position.x = curve.origin.x + std::remainder(collision_position.x - curve.origin.x,
                                              double(grid.geometry.nx));
    auto projection = curve.project(collision_position);
    auto fitted = curve.normal(projection.u);
    double alignment = dot(fitted, collision);
    if (alignment == 0) {
      use_curve = false;
      ++diagnostics.tangent_fallbacks;
    } else {
      orientation = alignment > 0 ? 1 : -1;
      diagnostics.reoriented_fits += orientation < 0;
      normal = fitted * orientation;
    }
  }
  double origin_distance = use_curve ? orientation * curve.project(curve.origin).distance : 0;

  // 一圈邻域只提供同一面连通分支内的候选，禁止跨气隙搬移库存。
  constexpr int radius = 1;
  std::vector<int> neighborhood;
  for (int dz = -radius; dz <= radius; ++dz)
    for (int dx = -radius; dx <= radius; ++dx)
      if (z + dz >= 0 && z + dz < grid.geometry.nz)
        neighborhood.push_back(grid.index(x + dx, z + dz));
  auto branch = update.connectivity.local_silicon_branch(original, neighborhood);

  // 原格夹角取零；其余候选沿用角度权重 0.3 与有符号距离评分。
  std::vector<RemovalCandidate> candidates{{original, -origin_distance, 0}};
  for (int dz = -radius; dz <= radius; ++dz)
    for (int dx = -radius; dx <= radius; ++dx) {
      int zz = z + dz;
      if (zz < 0 || zz >= grid.geometry.nz)
        continue;
      int cell = grid.index(x + dx, zz);
      if (cell == original || grid.material[cell] != Silicon || !grid.boundary[cell])
        continue;
      if (!branch.contains(cell)) {
        ++diagnostics.disconnected_donors;
        continue;
      }
      Vec offset{double(dx), double(dz)};
      double angle = std::acos(std::clamp(dot(unit(offset), normal), -1.0, 1.0));
      double distance = use_curve ? orientation * curve.project(curve.origin + offset).distance
                                  : dot(offset, collision);
      candidates.push_back({cell, normal_weight * angle - distance, dot(offset, offset)});
    }

  return update.apply(original, {std::move(branch), std::move(candidates)}, erase_cell);
}
}
