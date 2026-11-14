#pragma once
#include "vector.hpp"
#include "grid.hpp"
#include "common/random.hpp"

namespace etch::planar {
inline Distribution cosine_angle() {
    // The normalized 2-D flux density is cos(alpha) / 2.
    return [](Random& rng) { return std::asin(2 * uniform(rng) - 1); };
}
struct Particle { Vec position, direction; double energy; };
struct Source {
    Distribution angle, energy;
    Particle emit(Random& rng, const Geometry& geometry) const {
        double x = geometry.nx * uniform(rng);
        double alpha = angle(rng);
        double e = energy(rng);
        return {{x, 0.25}, {std::sin(alpha), std::cos(alpha)}, e};
    }
};
}

