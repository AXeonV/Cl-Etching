#pragma once
#include "grid.hpp"
#include "common/direction.hpp"
#include "common/particle.hpp"

namespace etch::spatial {
using Particle = ParticleState<Vec3>;

struct Source {
  DirectionDistribution direction;
  Distribution energy;

  // 源平面独立均匀采样 x、y，方向和能量由注入的分布函数生成。
  Particle emit(Random& rng, const Geometry& g) const {
    Vec3 position{g.nx * uniform(rng), g.ny * uniform(rng), 0.25};
    Vec3 d = direction(rng);
    return {position, d, energy(rng)};
  }
};
}
