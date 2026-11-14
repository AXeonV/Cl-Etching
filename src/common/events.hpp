#pragma once
namespace etch {
enum class RayEnd { Escaped, MissedSurface, PeriodicOrbit };
template<class Geometry> struct RunSettings {
    Geometry geometry;
    int fit_radius=3, max_bounces=0;
    bool reflections=true;
    unsigned long long seed=42;
};
struct Statistics {
    long long particles=0, silicon_hits=0, mask_hits=0;
    long long reflections=0, escaped=0, bounce_limit=0, deleted_cells=0;
    long long relocated_cells=0, pruned_si=0, pruned_cl=0;
};
}
