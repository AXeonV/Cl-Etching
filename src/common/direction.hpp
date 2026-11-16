#pragma once
#include "vector.hpp"
#include "random.hpp"
namespace etch {
using DirectionDistribution=std::function<Vec3(Random&)>;
inline DirectionDistribution normal_direction() { return [](Random&) { return Vec3{0,0,1}; }; }
inline DirectionDistribution cosine_direction() {
    return [](Random& rng) {
        double mu=std::sqrt(uniform(rng)),phi=2*std::acos(-1.0)*uniform(rng);
        double sine=std::sqrt(1-mu*mu);
        return Vec3{sine*std::cos(phi),sine*std::sin(phi),mu};
    };
}
// Gaussian intensity per solid angle: p(theta) is proportional to sin(theta)*I(theta).
inline DirectionDistribution gaussian_direction(double sigma) {
    return [sigma](Random& rng) {
        double theta;
        do {
            theta=sigma*std::sqrt(-2*std::log(1-uniform(rng)));
        } while(theta>=std::acos(-1.0)/2 || uniform(rng)>(theta==0?1:std::sin(theta)/theta));
        double phi=2*std::acos(-1.0)*uniform(rng);
        return Vec3{std::sin(theta)*std::cos(phi),std::sin(theta)*std::sin(phi),std::cos(theta)};
    };
}
}
