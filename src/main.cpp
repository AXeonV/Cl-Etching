#include "simulation.hpp"
#include <chrono>
#include <fstream>
#include <iostream>
#include <map>

int main(int argc, char** argv) {
    using namespace etch;
    std::map<std::string, std::string> args = {
        {"--model", "yield"}, {"--angular", "cosine"}, {"--angle", "gaussian"},
        {"--energy", "constant"}, {"--seed", "42"}, {"--output", "outputs/run"},
        {"--reflect", "1"}, {"--max-bounces", "0"},
        {"--nx", ""}, {"--nz", ""}, {"--surface-z", ""}, {"--opening", ""},
        {"--mask-thickness", ""}, {"--particles", ""}, {"--fwhm", ""},
        {"--energy-ev", ""}, {"--energy-sigma", ""}
    };
    for (int i = 1; i < argc; i += 2) args.at(argv[i]) = argv[i + 1];
    Settings settings{Geometry{std::stoi(args.at("--nx")), std::stoi(args.at("--nz")),
        std::stoi(args.at("--surface-z")), std::stoi(args.at("--mask-thickness")), std::stoi(args.at("--opening"))}};
    settings.seed = std::stoull(args.at("--seed"));
    settings.reflections = std::stoi(args.at("--reflect"));
    settings.max_bounces = std::stoi(args.at("--max-bounces"));
    const double pi = std::acos(-1.0);
    std::map<std::string, std::function<Distribution()>> angles = {
        {"gaussian", [&] { return gaussian(0, std::stod(args.at("--fwhm")) * pi / 180 / 2.354820045, -pi/2, pi/2); }},
        {"cosine", [] { return cosine_angle(); }}, {"normal", [] { return constant(0); }}
    };
    std::map<std::string, std::function<Distribution()>> energies = {
        {"constant", [&] { return constant(std::stod(args.at("--energy-ev"))); }},
        {"gaussian", [&] { return gaussian(std::stod(args.at("--energy-ev")),
            std::stod(args.at("--energy-sigma")), 0, std::numeric_limits<double>::infinity()); }}
    };
    std::map<std::string, AngularYield> laws = {{"cosine", cosine_yield}, {"enhanced", enhanced_yield}};
    std::map<std::string, std::function<std::unique_ptr<Reaction>()>> models = {
        {"yield", [&] { return std::make_unique<YieldReaction>(laws.at(args.at("--angular"))); }},
        {"chlorination", [&] { args["--angular"] = "not_used"; return std::make_unique<ChlorinationReaction>(); }}
    };
    Source source{angles.at(args.at("--angle"))(), energies.at(args.at("--energy"))()};
    Simulation simulation(settings, source, models.at(args.at("--model"))());
    auto start = std::chrono::steady_clock::now();
    long long particles = std::stoll(args.at("--particles"));
    long long batch = std::min(250000LL, particles);
    for (long long done = 0; done < particles; done += batch) {
        simulation.advance(std::min(batch, particles - done));
        std::cout << simulation.stats.particles << " particles, " << simulation.reaction->removed
                  << " Si removed, " << simulation.stats.deleted_cells << " cells" << std::endl;
    }
    double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    simulation.save(args.at("--output"), elapsed);
    std::ofstream config(args.at("--output") + "/arguments.txt");
    for (auto& [key, value] : args) config << key << " " << value << "\n";
    std::cout << "Saved " << args.at("--output") << " (" << elapsed << " s)" << std::endl;
}
