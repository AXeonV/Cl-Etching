#pragma once
#include "material.hpp"
#include "random.hpp"
#include <optional>

namespace etch {
// 将选中邻格的库存搬回原格，再把原格残量等概率分配到可用空位。
inline std::optional<CellRemoval> transfer_inventory(std::vector<Sites>& sites, int original,
  int selected, const std::vector<int>& recipients, int capacity, Random& rng) {
  Sites pending = sites[original];
  int residual = site_count(pending), available = 0;
  // 先按搬移后的布局核算容量；原格将接收选中邻格的库存。
  for (int cell : recipients) {
    int source = cell == original ? selected : cell;
    available += capacity - site_count(sites[source]);
  }
  // 容量不足时不修改库存，也不消耗随机数；由上层延后几何更新。
  if (available < residual)
    return std::nullopt;

  if (selected != original)
    sites[original] = sites[selected];
  sites[selected] = {};
  CellRemoval removal{selected, selected != original};
  removal.redistributed_si = residual;
  for (int state = 1; state < 4; ++state)
    removal.redistributed_cl += state * pending[state];

  std::vector<int> free_slots;
  for (int cell : recipients)
    free_slots.push_back(capacity - site_count(sites[cell]));
  while (residual > 0) {
    // 对空位等概率抽样，接收格的概率正比于剩余容量。
    int target = static_cast<int>(uniform(rng) * available), destination = 0;
    while (target >= free_slots[destination])
      target -= free_slots[destination++];

    // 再按残余各组分的数量抽取状态，保持 Si 及已结合 Cl 的总数。
    target = static_cast<int>(uniform(rng) * residual);
    int state = 0;
    while (target >= pending[state])
      target -= pending[state++];

    --pending[state];
    ++sites[recipients[destination]][state];
    --free_slots[destination];
    --residual;
    --available;
  }
  return removal;
}
}
