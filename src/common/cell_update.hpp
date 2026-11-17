#pragma once
#include "cell_inventory.hpp"
#include "connectivity.hpp"
#include "events.hpp"

namespace etch {
struct CellNeighborhood {
    SiliconBranch branch;
    std::vector<RemovalCandidate> candidates;
};

// Dimensional backends supply geometry scores and an erase/cache callback.
template<class Grid> class CellUpdate {
    Grid& grid;
    CellUpdateDiagnostics& diagnostics;
    Random rng;
public:
    RemovalConnectivity<Grid> connectivity;

    CellUpdate(Grid& grid,CellUpdateDiagnostics& diagnostics,unsigned long long seed)
        :grid(grid),diagnostics(diagnostics),rng(seed),connectivity(grid) {}

    bool isolated(int original) const {
        if(grid.is_anchor(original))return false;
        for(int neighbor:grid.face_neighbors(original))
            if(neighbor>=0&&neighbor!=original&&grid.material[neighbor]!=Gas)return false;
        return true;
    }

    template<class Erase> CellRemoval detach(int cell,Erase erase) {
        CellRemoval removal{cell,false};removal.reaction_cells=0;
        release_inventory(grid.sites,cell,removal);
        erase(cell);
        return removal;
    }

    template<class Erase> CellRemoval apply(int original,CellNeighborhood neighborhood,Erase erase) {
        auto choice=connectivity.select(neighborhood.candidates,rng);
        int selected=choice.cell;
        std::vector<int> recipients;
        for(int cell:neighborhood.branch.cells)if(cell!=selected)recipients.push_back(cell);
        auto transferred=transfer_inventory(grid.sites,original,selected,recipients,grid.geometry.atoms_per_cell,rng);
        if(!transferred) {
            ++diagnostics.deferred_updates;
            CellRemoval pending{original,false};pending.reaction_cells=0;
            return pending;
        }
        auto removal=*transferred;
        diagnostics.redirected_removals+=choice.redirected;
        diagnostics.unavoidable_splits+=!choice.connected;
        bool anchor_removed=grid.is_anchor(selected);
        erase(selected);
        if(!choice.connected||anchor_removed)for(int cell:connectivity.detached_components(selected)) {
            release_inventory(grid.sites,cell,removal);
            erase(cell);
        }
        return removal;
    }
};
}
