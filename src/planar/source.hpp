#pragma once
#include "grid.hpp"
#include "common/direction.hpp"
#include "common/particle.hpp"

namespace etch::planar {
using Particle = ParticleState<Vec>;

struct Source {
  DirectionDistribution direction;
  Distribution energy;

  // 源线均匀采样位置，保留完整三分量方向供碰撞与反射使用。
  Particle emit(Random& rng, const Geometry& geometry) const {
    double x = geometry.nx * uniform(rng);
    Vec3 d = direction(rng);
    return {{x, 0.25}, d, energy(rng)};
  }
};
}
