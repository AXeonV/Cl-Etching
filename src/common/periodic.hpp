#pragma once
#include <cmath>
#include <algorithm>
#include <cstdint>
#include <numeric>
#include <limits>

namespace etch {
// 水平射线在两个横向位移都为整数周期时闭合。
// 约分浮点方向的精确二进制有效位，直接计算共同周期。
inline double horizontal_ray_period(double dx, double dy, int nx, int ny) {
  if (dx == 0)
    return ny / std::abs(dy);
  if (dy == 0)
    return nx / std::abs(dx);
  auto parts = [](double value, int& exponent) {
    double mantissa = std::frexp(std::abs(value), &exponent);
    uint64_t significant = static_cast<uint64_t>(std::ldexp(mantissa, 53));
    exponent -= 53;
    while ((significant & 1) == 0) {
      significant >>= 1;
      ++exponent;
    }
    return significant;
  };
  int ex, ey;
  uint64_t sx = parts(dx, ex), sy = parts(dy, ey), x = nx, y = ny;
  auto cancel = [](uint64_t& a, uint64_t& b) {
    uint64_t g = std::gcd(a, b);
    a /= g;
    b /= g;
  };

  // 分别约分有效位与网格周期，再补回二进制指数差。
  cancel(sx, sy);
  cancel(x, y);
  cancel(sx, x);
  cancel(sy, y);
  int shift = ex - ey;
  while (shift > 0 && (x & 1) == 0) {
    x >>= 1;
    --shift;
  }
  while (shift < 0 && (y & 1) == 0) {
    y >>= 1;
    ++shift;
  }
  long double cycles = std::ldexp(static_cast<long double>(sx) * y, std::max(shift, 0));
  return static_cast<double>(cycles * nx / std::abs(dx));
}
}
