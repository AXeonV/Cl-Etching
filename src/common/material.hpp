#pragma once
#include <array>
#include <cstdint>
#include <numeric>
#include <vector>

namespace etch {
enum Material : uint8_t { Gas, Silicon, Mask };
using Sites = std::array<uint16_t,4>;
struct CellRemoval {
    int cell;bool relocated;
    long long detached_cells=0,detached_si=0,detached_cl=0;
    long long redistributed_si=0,redistributed_cl=0;
    int reaction_cells=1;
};
inline void release_inventory(std::vector<Sites>& sites,int cell,CellRemoval& removal) {
    for(int state=0;state<4;++state) {
        removal.detached_si+=sites[cell][state];removal.detached_cl+=state*sites[cell][state];
    }
    ++removal.detached_cells;sites[cell]={};
}
inline int site_count(const Sites& sites) { return std::accumulate(sites.begin(),sites.end(),0); }
inline long long silicon_inventory(const std::vector<Sites>& sites) {
    long long count=0;for(const auto& cell:sites)count+=site_count(cell);return count;
}
inline long long chlorine_inventory(const std::vector<Sites>& sites) {
    long long count=0;for(const auto& cell:sites)for(int state=1;state<4;++state)count+=state*cell[state];return count;
}
}
