#pragma once
#include <array>
#include <cstdint>
#include <numeric>
#include <vector>

namespace etch {
enum Material : uint8_t { Gas, Silicon, Mask };
using Sites = std::array<uint16_t,4>;
inline int site_count(const Sites& sites) { return std::accumulate(sites.begin(),sites.end(),0); }
inline long long silicon_inventory(const std::vector<Sites>& sites) {
    long long count=0;for(const auto& cell:sites)count+=site_count(cell);return count;
}
inline long long chlorine_inventory(const std::vector<Sites>& sites) {
    long long count=0;for(const auto& cell:sites)for(int state=1;state<4;++state)count+=state*cell[state];return count;
}
}
