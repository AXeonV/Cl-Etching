#pragma once
#include "cell_inventory.hpp"
#include "connectivity.hpp"
#include "events.hpp"

namespace etch {
struct CellNeighborhood {
  SiliconBranch branch;
  std::vector<RemovalCandidate> candidates;
};

// 后端提供几何评分和删格回调；公共流程负责守恒迁移、碎片释放与统计。
template <class Grid> class CellUpdate {
  Grid& grid;
  CellUpdateDiagnostics& diagnostics;
  Random rng;

public:
  RemovalConnectivity<Grid> connectivity;

  CellUpdate(Grid& grid, CellUpdateDiagnostics& diagnostics, unsigned long long seed)
      : grid(grid), diagnostics(diagnostics), rng(seed), connectivity(grid) {}

  bool isolated(int original) const {
    if (grid.is_anchor(original))
      return false;
    for (int neighbor : grid.face_neighbors(original))
      if (neighbor >= 0 && neighbor != original && grid.material[neighbor] != Gas)
        return false;
    return true;
  }

  template <class Erase> CellRemoval detach(int cell, Erase erase) {
    CellRemoval removal{cell, false};
    removal.reaction_cells = 0;
    release_inventory(grid.sites, cell, removal);
    erase(cell);
    return removal;
  }

  template <class Erase>
  CellRemoval apply(int original, CellNeighborhood neighborhood, Erase erase) {
    auto choice = connectivity.select(neighborhood.candidates, rng);
    int selected = choice.cell;

    // 删除格不接收残量；其余接收格均来自原格的局部面连通分支。
    std::vector<int> recipients;
    for (int cell : neighborhood.branch.cells)
      if (cell != selected)
        recipients.push_back(cell);

    auto transferred = transfer_inventory(
      grid.sites, original, selected, recipients, grid.geometry.atoms_per_cell, rng);
    if (!transferred) {
      ++diagnostics.deferred_updates;
      CellRemoval pending{original, false};
      pending.reaction_cells = 0;
      return pending;
    }

    auto removal = *transferred;
    diagnostics.redirected_removals += choice.redirected;
    diagnostics.unavoidable_splits += !choice.connected;

    // 只有库存搬移成功后才改材料标签，并通过后端回调使几何缓存失效。
    bool anchor_removed = grid.is_anchor(selected);
    erase(selected);

    // 失去锚点连接的分量作为固体碎片释放，不计入化学产物。
    if (!choice.connected || anchor_removed)
      for (int cell : connectivity.detached_components(selected)) {
        release_inventory(grid.sites, cell, removal);
        erase(cell);
      }
    return removal;
  }
};
}
