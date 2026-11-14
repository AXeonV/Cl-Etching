#include "reaction.hpp"
#include <algorithm>

namespace etch {
void YieldReaction::apply(Sites& sites, double mu, double capture_probability, Random& rng) {
    // F is the yield per incident ion; this event is conditional on capture.
    double expected = angular(mu) / capture_probability;
    int whole = static_cast<int>(expected);
    int count = whole + (uniform(rng) < expected - whole);
    count = std::min(count, site_count(sites));
    sites[0] -= count;
    ++captures;
    removed += count;
}
void ChlorinationReaction::apply(Sites& sites, double, double, Random& rng) {
    int target = static_cast<int>(uniform(rng) * site_count(sites));
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
