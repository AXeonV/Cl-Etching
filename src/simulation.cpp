#include "simulation.hpp"
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>

namespace etch {
Simulation::Simulation(Settings s, Source src, std::unique_ptr<Reaction> rx, CaptureProbability p)
    : settings(s), grid(s.geometry), surface(grid, s.fit_radius), transport(grid, surface),
      source(std::move(src)), capture(std::move(p)), reaction(std::move(rx)),
      source_rng(s.seed), event_rng(s.seed + 1000003), primary_hits(grid.material.size()),
      secondary_hits(grid.material.size()), initial_si(grid.silicon_inventory()) {}

void Simulation::advance(long long count) {
    for (long long n = 0; n < count; ++n) {
        Particle particle = source.emit(source_rng, grid.geometry.nx);
        ++stats.particles;
        for (int bounce = 0; ; ++bounce) {
            auto hit = transport.next_hit(particle);
            if (!hit) {
                if (transport.last_end == RayEnd::Escaped) ++stats.escaped;
                break;
            }
            int cell = grid.index(hit->x, hit->z);
            double mu = -dot(particle.direction, hit->normal);
            bool reflect = false;
            if (grid.material[cell] == Mask) {
                ++stats.mask_hits;
                reflect = true;
            } else {
                ++stats.silicon_hits;
                (bounce == 0 ? primary_hits : secondary_hits)[cell]++;
                if (uniform(event_rng) < capture(mu)) {
                    reaction->apply(grid, cell, mu, event_rng);
                    if (grid.remaining(cell) < grid.geometry.deletion_threshold) {
                        stats.pruned_si += grid.remaining(cell);
                        for (int state = 1; state < 4; ++state) stats.pruned_cl += state * grid.sites[cell][state];
                        grid.sites[cell] = {};
                        grid.remove_empty_cell(hit->x, hit->z);
                        surface.invalidate(hit->x, hit->z);
                        ++stats.deleted_cells;
                    }
                    break;
                }
                reflect = true;
            }
            if (!reflect || !settings.reflections) break;
            if (settings.max_bounces > 0 && bounce == settings.max_bounces) { ++stats.bounce_limit; break; }
            particle.position = hit->position + hit->normal * 1e-5;
            particle.position.x -= grid.geometry.nx * std::floor(particle.position.x / grid.geometry.nx);
            particle.direction = reflected(particle.direction, hit->normal);
            ++stats.reflections;
        }
    }
}

template<class T> void write_binary(const std::string& path, const std::vector<T>& data) {
    std::ofstream file(path, std::ios::binary);
    file.write(reinterpret_cast<const char*>(data.data()), data.size() * sizeof(T));
}
void Simulation::save(const std::string& directory, double elapsed_seconds) {
    std::filesystem::create_directories(directory);
    write_binary(directory + "/material.bin", grid.material);
    write_binary(directory + "/sites.bin", grid.sites);
    write_binary(directory + "/primary_hits.bin", primary_hits);
    write_binary(directory + "/secondary_hits.bin", secondary_hits);
    std::ofstream curves(directory + "/surface.csv");
    curves << "x,z,tx,tz,nx,nz,a,b,c,cell\n" << std::setprecision(15);
    for (int cell : grid.surface_cells) {
        const auto& c = surface.curve(cell % grid.geometry.nx, cell / grid.geometry.nx);
        curves << c.origin.x << ',' << c.origin.z << ',' << c.tangent.x << ',' << c.tangent.z << ','
               << c.outward.x << ',' << c.outward.z << ',' << c.a << ',' << c.b << ',' << c.c << ',' << cell << '\n';
    }
    auto g = grid.geometry;
    std::ofstream f(directory + "/summary.json");
    f << std::setprecision(12)
      << "{\n  \"model\": \"" << reaction->name() << "\",\n"
      << "  \"nx\": " << g.nx << ", \"nz\": " << g.nz << ",\n"
      << "  \"cell_nm\": " << g.cell_nm << ", \"surface_z\": " << g.surface_z << ",\n"
      << "  \"opening_cells\": " << g.opening << ", \"mask_cells\": " << g.mask_thickness << ",\n"
      << "  \"atoms_per_cell\": " << g.atoms_per_cell << ", \"seed\": " << settings.seed << ",\n"
      << "  \"deletion_threshold\": " << g.deletion_threshold << ",\n"
      << "  \"fit_radius\": " << settings.fit_radius << ", \"max_bounces\": " << settings.max_bounces << ",\n"
      << "  \"reflections_enabled\": " << (settings.reflections ? "true" : "false") << ",\n"
      << "  \"particles\": " << stats.particles << ", \"silicon_hits\": " << stats.silicon_hits << ",\n"
      << "  \"mask_hits\": " << stats.mask_hits << ", \"reflections\": " << stats.reflections << ",\n"
      << "  \"escaped\": " << stats.escaped << ", \"bounce_limit\": " << stats.bounce_limit << ",\n"
      << "  \"interior_exits\": " << transport.interior_exits << ",\n"
      << "  \"normal_mapped_hits\": " << transport.normal_mapped_hits << ",\n"
      << "  \"mapping_failures\": " << transport.mapping_failures << ",\n"
      << "  \"no_local_intersections\": " << transport.no_local_intersections << ",\n"
      << "  \"pruned_si\": " << stats.pruned_si << ", \"pruned_cl\": " << stats.pruned_cl << ",\n"
      << "  \"periodic_orbits\": " << transport.periodic_orbits << ",\n"
      << "  \"rank_deficient_fits\": " << surface.rank_deficient_fits << ",\n"
      << "  \"captures\": " << reaction->captures << ", \"removed_si\": " << reaction->removed << ",\n"
      << "  \"sicl4_products\": " << reaction->products_sicl4 << ",\n"
      << "  \"initial_si\": " << initial_si << ", \"remaining_si\": " << grid.silicon_inventory() << ",\n"
      << "  \"surface_cl\": " << grid.chlorine_inventory() << ",\n"
      << "  \"deleted_cells\": " << stats.deleted_cells << ", \"fits\": " << surface.fit_count << ",\n"
      << "  \"elapsed_seconds\": " << elapsed_seconds << "\n}\n";
}
}

