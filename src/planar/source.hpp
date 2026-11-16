#pragma once
#include "grid.hpp"
#include "common/direction.hpp"
#include "common/particle.hpp"

namespace etch::planar {
using Particle=ParticleState<Vec>;
struct Source {
    DirectionDistribution direction;
    Distribution energy;
    Particle emit(Random& rng, const Geometry& geometry) const {
        double x = geometry.nx * uniform(rng);
        Vec3 d = direction(rng);
        return {{x, 0.25}, d, energy(rng)};
    }
};
}
