#pragma once
#include "vector.hpp"
namespace etch {
template<class Position> struct ParticleState {
    Position position;
    Vec3 direction;
    double energy;
};
template<class Position> struct SurfaceHit {
    Position position;
    Vec3 normal;
    int cell;
};
}
