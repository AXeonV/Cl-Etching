#pragma once
#include "material.hpp"
#include "random.hpp"
#include <algorithm>
#include <cstdint>
#include <vector>

namespace etch {
struct RemovalCandidate {
  int cell;
  double score;
  double distance_squared = 0;
};

struct RemovalChoice {
  int cell;
  bool connected, redirected;
};

struct SiliconBranch {
  std::vector<int> cells;

  bool contains(int cell) const {
    return std::find(cells.begin(), cells.end(), cell) != cells.end();
  }
};

// 网格提供带横向周期边界的面邻接，-1 表示超出开放边界。
template <class Grid> class RemovalConnectivity {
  const Grid& grid;
  // 用搜索代次复用访问标记，避免每次判定都清空整个网格。
  std::vector<uint64_t> visited;
  std::vector<uint8_t> owner;
  std::array<std::vector<int>, 6> fronts;
  uint64_t generation = 0;
  std::vector<int> queue;

public:
  explicit RemovalConnectivity(const Grid& g)
      : grid(g), visited(g.material.size()), owner(g.material.size()) {}

  SiliconBranch local_silicon_branch(int original, const std::vector<int>& neighborhood) const {
    SiliconBranch branch{{original}};
    for (size_t i = 0; i < branch.cells.size(); ++i)
      for (int neighbor : grid.face_neighbors(branch.cells[i])) {
        if (neighbor < 0 || grid.material[neighbor] != Silicon || branch.contains(neighbor))
          continue;
        if (std::find(neighborhood.begin(), neighborhood.end(), neighbor) == neighborhood.end())
          continue;
        branch.cells.push_back(neighbor);
      }
    return branch;
  }

  // 虚拟删除后，其余固体面邻格仍互相可达，才保持原有固体连接。
  // 检查可超出局部候选范围，避免误判远处绕行或跨周期边界的连接。
  bool preserves_connection(int removed) {
    std::vector<int> targets;
    for (int neighbor : grid.face_neighbors(removed)) {
      if (neighbor < 0 || neighbor == removed || grid.material[neighbor] == Gas)
        continue;
      if (std::find(targets.begin(), targets.end(), neighbor) == targets.end())
        targets.push_back(neighbor);
    }
    if (targets.size() < 2)
      return true;

    // 从各个存活面邻格同时扩展，波前相遇时合并连通分组。
    ++generation;
    std::array<int, 6> parents{};
    std::array<size_t, 6> heads{};
    int count = static_cast<int>(targets.size()), groups = count;
    for (int i = 0; i < count; ++i) {
      parents[i] = i;
      fronts[i].clear();
      fronts[i].push_back(targets[i]);
      visited[targets[i]] = generation;
      owner[targets[i]] = static_cast<uint8_t>(i);
    }

    auto root = [&](int id) {
      while (parents[id] != id) {
        parents[id] = parents[parents[id]];
        id = parents[id];
      }
      return id;
    };

    for (;;) {
      for (int i = 0; i < count; ++i)
        if (heads[i] < fronts[i].size()) {
          int cell = fronts[i][heads[i]++];
          for (int neighbor : grid.face_neighbors(cell)) {
            if (neighbor < 0 || neighbor == removed || grid.material[neighbor] == Gas)
              continue;
            if (visited[neighbor] != generation) {
              visited[neighbor] = generation;
              owner[neighbor] = static_cast<uint8_t>(i);
              fronts[i].push_back(neighbor);
            } else {
              int a = root(i), b = root(owner[neighbor]);
              if (a != b) {
                parents[a] = b;
                if (--groups == 1)
                  return true;
              }
            }
          }
        }

      // 某分组的全部波前耗尽而仍未合并，说明它已与其余分组断开。
      std::array<bool, 6> active{};
      for (int i = 0; i < count; ++i)
        if (heads[i] < fronts[i].size())
          active[root(i)] = true;
      for (int i = 0; i < count; ++i)
        if (parents[i] == i && !active[i])
          return false;
    }
  }

  // 掩膜或底部硅提供锚点，搜到任意锚点即可保留整个固体分量。
  std::vector<int> detached_components(int removed) {
    std::vector<int> detached;
    for (int seed : grid.face_neighbors(removed)) {
      if (seed < 0 || grid.material[seed] == Gas ||
          std::find(detached.begin(), detached.end(), seed) != detached.end())
        continue;
      ++generation;
      queue.clear();
      queue.push_back(seed);
      visited[seed] = generation;
      bool anchored = false;
      for (size_t i = 0; i < queue.size() && !anchored; ++i) {
        int cell = queue[i];
        if (grid.is_anchor(cell)) {
          anchored = true;
          break;
        }
        for (int neighbor : grid.face_neighbors(cell)) {
          if (neighbor < 0 || grid.material[neighbor] == Gas || visited[neighbor] == generation)
            continue;
          visited[neighbor] = generation;
          queue.push_back(neighbor);
        }
      }
      if (!anchored)
        detached.insert(detached.end(), queue.begin(), queue.end());
    }
    return detached;
  }

  RemovalChoice select(std::vector<RemovalCandidate>& candidates, Random& rng) {
    auto rank = [](const auto& a, const auto& b) {
      return a.score < b.score || (a.score == b.score && a.distance_squared < b.distance_squared);
    };

    // 周期映射可能产生同格候选；先去重，再按几何评分和位移排序。
    std::sort(candidates.begin(), candidates.end(), [&](const auto& a, const auto& b) {
      return a.cell < b.cell || (a.cell == b.cell && rank(a, b));
    });
    candidates.erase(std::unique(candidates.begin(), candidates.end(),
                       [](auto a, auto b) { return a.cell == b.cell; }),
      candidates.end());
    std::sort(candidates.begin(), candidates.end(), rank);

    size_t i = 0;
    while (i < candidates.size()) {
      size_t end = i + 1;
      while (end < candidates.size() && !rank(candidates[i], candidates[end]))
        ++end;
      // 完全等价的安全候选随机选择，避免固定的格点方向偏置。
      std::vector<int> safe;
      for (size_t j = i; j < end; ++j)
        if (preserves_connection(candidates[j].cell))
          safe.push_back(candidates[j].cell);
      if (!safe.empty()) {
        int selected =
          safe.size() == 1 ? safe.front() : safe[static_cast<int>(uniform(rng) * safe.size())];
        return {selected, true, i > 0};
      }
      i = end;
    }

    // 若所有候选都会断裂，仍按最佳评分推进，交由碎片路径处理。
    size_t end = 1;
    while (end < candidates.size() && !rank(candidates.front(), candidates[end]))
      ++end;
    size_t selected = end == 1 ? 0 : static_cast<size_t>(uniform(rng) * end);
    return {candidates[selected].cell, false, false};
  }
};
}
