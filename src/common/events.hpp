#pragma once
#include <array>
namespace etch {
enum class RayEnd { Escaped, MissedSurface, PeriodicOrbit };
struct TrackingDiagnostics {
    RayEnd last_end=RayEnd::Escaped;
    long long interior_exits=0,periodic_orbits=0;
    long long normal_mapped_hits=0,mapping_failures=0,off_normal_mapped_hits=0;
    double max_mapping_offset=0;
};
template<class Geometry> struct RunSettings {
    Geometry geometry;
    int fit_radius=3, max_bounces=0, surface_subdivisions=1;
    bool reflections=true;
    unsigned long long seed=42;
    double smoothing_weight=0.3;
};
struct CellUpdateDiagnostics {
    long long undefined_fit_removals=0,disconnected_donors=0;
    long long redirected_removals=0,unavoidable_splits=0;
    long long deferred_updates=0;
    long long reoriented_fits=0,tangent_fallbacks=0;
};
struct AngleBin { long long primary=0,secondary=0,captures=0,removed=0; };
struct Statistics {
    std::array<AngleBin,90> angles{};
    long long particles=0, silicon_hits=0, mask_hits=0;
    long long reflections=0, escaped=0, bounce_limit=0, deleted_cells=0;
    long long relocated_cells=0;
    long long detached_cells=0,detached_si=0,detached_cl=0;
    long long redistributed_si=0,redistributed_cl=0;
};
}
