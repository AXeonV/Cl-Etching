#include "reaction.hpp"
#include <algorithm>

namespace etch {
void YieldReaction::apply(Grid& grid, int cell, double mu, Random& rng) {
    // The paper computes the etching yield after the particle enters the surface layer.
    double expected = angular(mu);
    int whole = static_cast<int>(expected);
    int count = whole + (uniform(rng) < expected - whole);
    count = std::min(count, grid.remaining(cell));
    grid.sites[cell][0] -= count;
    ++captures;
    removed += count;
}
void ChlorinationReaction::apply(Grid& grid, int cell, double, Random& rng) {
    auto& sites = grid.sites[cell];
    int target = static_cast<int>(uniform(rng) * grid.remaining(cell));
    int state = 0;
    while (target >= sites[state]) target -= sites[state++];
    --sites[state];
    ++captures;
    if (state == 3) {
        ++removed;
        ++products_sicl4;
    } else ++sites[state + 1];
}
}
