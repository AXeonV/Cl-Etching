#pragma once
#include <random>
#include <functional>
#include <cmath>

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
    do {
      value = distribution(rng);
    } while (value < lower || value > upper);
    return value;
  };
}
}
