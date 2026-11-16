#pragma once
#include "grid.hpp"
#include "common/direction.hpp"
#include "common/particle.hpp"
namespace etch::spatial {
using Particle=ParticleState<Vec3>;
struct Source {
    DirectionDistribution direction;
    Distribution energy;
    Particle emit(Random& rng,const Geometry& g) const {
        Vec3 position{g.nx*uniform(rng),g.ny*uniform(rng),0.25};
        Vec3 d=direction(rng);return {position,d,energy(rng)};
    }
};
}
