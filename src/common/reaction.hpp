#pragma once
#include "material.hpp"
#include "random.hpp"
#include <functional>
#include <string>

namespace etch {
using AngularYield = std::function<double(double)>;
using CaptureProbability = std::function<double(double)>;

inline double cosine_yield(double mu) {
  return mu;
}

inline double enhanced_yield(double mu) {
  return (5 - 4 * mu * mu) * mu;
}

inline double cosine_capture(double mu) {
  return mu;
}

// 仅处理已捕获粒子的反应；几何删除和固体碎片账目由更新器负责。
class Reaction {
public:
  virtual ~Reaction() = default;
  virtual void apply(Sites& sites, double mu, double capture_probability, Random& rng) = 0;
  virtual std::string name() const = 0;
  long long captures = 0, removed = 0, products_sicl4 = 0;
};

class YieldReaction : public Reaction {
  AngularYield angular;

public:
  explicit YieldReaction(AngularYield angular = cosine_yield) : angular(angular) {}

  void apply(Sites& sites, double mu, double capture_probability, Random& rng) override;

  std::string name() const override {
    return "yield";
  }
};

class ChlorinationReaction : public Reaction {
public:
  void apply(Sites& sites, double mu, double capture_probability, Random& rng) override;

  std::string name() const override {
    return "chlorination";
  }
};
}
