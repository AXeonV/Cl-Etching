#pragma once
#include "reaction.hpp"
#include "cell_update.hpp"
#include "transport.hpp"
#include <memory>
#include <string>

namespace etch {
struct Settings {
    Geometry geometry;
    int fit_radius = 3, max_bounces = 0;
    bool reflections = true;
    unsigned long long seed = 42;
};
struct Statistics {
    long long particles = 0, silicon_hits = 0, mask_hits = 0;
    long long reflections = 0, escaped = 0, bounce_limit = 0, deleted_cells = 0;
    long long relocated_cells = 0;
    long long pruned_si = 0, pruned_cl = 0;
};
class Simulation {
public:
    Settings settings;
    Grid grid;
    Surface surface;
    Transport transport;
    CellUpdater cell_updater;
    Source source;
    CaptureProbability capture;
    std::unique_ptr<Reaction> reaction;
    Random source_rng, event_rng;
    Statistics stats;
    std::vector<uint32_t> primary_hits, secondary_hits;
    long long initial_si;
    Simulation(Settings settings, Source source, std::unique_ptr<Reaction> reaction,
               CaptureProbability capture = cosine_capture);
    void advance(long long count);
    void save(const std::string& directory, double elapsed_seconds);
};
}

