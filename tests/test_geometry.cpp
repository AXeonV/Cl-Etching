#include "simulation.hpp"
#include <iostream>

using namespace etch;

int main() {
    // Software diagnostic fixture, independent of the paper's figure inputs.
    Geometry geometry{256, 512, 60, 24, 72};
    Grid flat(geometry); Surface surface(flat); Transport transport(flat, surface);
    int wrong_cells = 0;
    for (int i = 0; i < 100; ++i) {
        double x = 125 + (i + 0.5) / 100;
        auto hit = transport.next_hit({{x, 1}, {0, 1}, 100});
        if (!hit || hit->x != static_cast<int>(std::floor(x))) ++wrong_cells;
    }
    std::cout << (wrong_cells ? "FAIL " : "PASS ") << "flat-surface hit ownership: "
              << wrong_cells << "/100 assigned to a different cell\n";

    Settings settings{geometry};
    settings.max_bounces = 500;
    double sigma = 10 * std::acos(-1.0) / 180 / 2.354820045;
    Source source{gaussian(0, sigma, -std::acos(-1.0) / 2, std::acos(-1.0) / 2), constant(100)};
    Simulation simulation(settings, source, std::make_unique<YieldReaction>());
    simulation.advance(100000);
    std::cout << (simulation.transport.interior_exits ? "FAIL " : "PASS ") << "surface interception: "
              << simulation.transport.interior_exits << "/100000 rays exited through the substrate\n";
    return wrong_cells || simulation.transport.interior_exits ? 1 : 0;
}
