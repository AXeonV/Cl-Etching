#pragma once
#include "vector.hpp"

namespace etch {
// 核心位置以格长为单位；导出时通过 cell_nm 转换为 nm。
template <class Position> struct ParticleState {
  Position position;
  // 二维与三维均保存三分量单位方向，二维只投影位置推进。
  Vec3 direction;
  double energy;
};

template <class Position> struct SurfaceHit {
  Position position;
  Vec3 normal;
  int cell;
};
}
