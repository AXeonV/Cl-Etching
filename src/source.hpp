#pragma once
#include "vector.hpp"
#include <functional>
#include <random>

namespace etch {
using Random = std::mt19937_64;
using Distribution = std::function<double(Random&)>;
inline double uniform(Random& rng) {
    return std::generate_canonical<double, 53>(rng);
}
inline Distribution constant(double value) {
    return [value](Random&) { return value; };
}
inline Distribution gaussian(double mean, double sigma, double lower, double upper) {
    return [=](Random& rng) {
        std::normal_distribution<double> distribution(mean, sigma);
        double value;
        do { value = distribution(rng); } while (value < lower || value > upper);
        return value;
    };
}
inline Distribution cosine_angle() {
    // The normalized 2-D flux density is cos(alpha) / 2.
    return [](Random& rng) { return std::asin(2 * uniform(rng) - 1); };
}
struct Particle { Vec position, direction; double energy; };
struct Source {
    Distribution angle, energy;
    Particle emit(Random& rng, double width) const {
        double x = width * uniform(rng);
        double alpha = angle(rng);
        double e = energy(rng);
        return {{x, 0.25}, {std::sin(alpha), std::cos(alpha)}, e};
    }
};
}

