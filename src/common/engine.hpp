#pragma once
#include "reaction.hpp"
#include "events.hpp"
#include "output.hpp"
#include <filesystem>
#include <memory>
#include <algorithm>

namespace etch {
template<class Backend> class Engine {
public:
    using Settings=RunSettings<typename Backend::Geometry>;
    Settings settings;
    typename Backend::Grid grid;
    typename Backend::Surface surface;
    typename Backend::Transport transport;
    typename Backend::CellUpdater cell_updater;
    typename Backend::Source source;
    CaptureProbability capture;
    std::unique_ptr<Reaction> reaction;
    Random source_rng,event_rng;
    Statistics stats;
    std::vector<uint32_t> primary_hits,secondary_hits;
    long long initial_si;

    Engine(Settings s,typename Backend::Source src,std::unique_ptr<Reaction> rx,
           CaptureProbability p=cosine_capture)
        :settings(s),grid(s.geometry),surface(grid,s.fit_radius,s.surface_subdivisions),transport(grid,surface),
         cell_updater(grid,surface),source(std::move(src)),capture(std::move(p)),reaction(std::move(rx)),
         source_rng(s.seed),event_rng(s.seed+1000003),primary_hits(grid.material.size()),
         secondary_hits(grid.material.size()),initial_si(grid.silicon_inventory()) {}

    void advance(long long count) {
        for(long long n=0;n<count;++n) {
            auto particle=source.emit(source_rng,grid.geometry);
            ++stats.particles;
            for(int bounce=0;;++bounce) {
                auto hit=transport.next_hit(particle);
                if(!hit) { if(transport.last_end==RayEnd::Escaped)++stats.escaped;break; }
                int cell=hit->cell;
                double mu=-dot(particle.direction,hit->normal);
                if(grid.material[cell]==Mask) ++stats.mask_hits;
                else {
                    ++stats.silicon_hits;
                    int angle=std::min(89,int(std::acos(std::clamp(mu,0.0,1.0))*180/std::acos(-1.0)));
                    auto& bin=stats.angles[angle];++(bounce==0?bin.primary:bin.secondary);
                    (bounce==0?primary_hits:secondary_hits)[cell]++;
                    double probability=capture(mu);
                    if(uniform(event_rng)<probability) {
                        long long before=reaction->removed;
                        reaction->apply(grid.sites[cell],mu,probability,event_rng);
                        ++bin.captures;bin.removed+=reaction->removed-before;
                        if(grid.remaining(cell)<grid.geometry.deletion_threshold) {
                            auto removal=cell_updater.remove(cell);
                            for(int state=0;state<4;++state) {
                                stats.pruned_si+=removal.residual[state];
                                stats.pruned_cl+=state*removal.residual[state];
                            }
                            stats.relocated_cells+=removal.relocated;
                            ++stats.deleted_cells;
                        }
                        break;
                    }
                }
                if(!settings.reflections)break;
                if(settings.max_bounces>0 && bounce==settings.max_bounces) { ++stats.bounce_limit;break; }
                particle.position=hit->position;
                grid.wrap_position(particle.position);
                particle.direction=reflected(particle.direction,hit->normal);
                ++stats.reflections;
            }
        }
    }

    void save(const std::string& directory,double elapsed_seconds) {
        std::filesystem::create_directories(directory);
        write_binary(directory+"/material.bin",grid.material);
        write_binary(directory+"/sites.bin",grid.sites);
        write_binary(directory+"/primary_hits.bin",primary_hits);
        write_binary(directory+"/secondary_hits.bin",secondary_hits);
        Backend::write_surface(surface,grid,directory+"/surface.csv");
        Backend::write_field(surface,grid,directory+"/surface_nodes.bin");
        std::ofstream angles(directory+"/collision_angles.csv");
        angles<<"angle_lower_deg,angle_upper_deg,primary_hits,secondary_hits,captures,removed_si\n";
        for(size_t i=0;i<stats.angles.size();++i) {
            const auto& bin=stats.angles[i];angles<<i<<','<<i+1<<','<<bin.primary<<','<<bin.secondary<<','<<bin.captures<<','<<bin.removed<<'\n';
        }
        Summary out(directory+"/summary.json");
        out.field("model",reaction->name());
        out.field("surface_representation",std::string("hermite"));
        Backend::write_geometry(out,grid.geometry);
        out.field("atoms_per_cell",grid.geometry.atoms_per_cell);
        out.field("seed",settings.seed);out.field("deletion_threshold",grid.geometry.deletion_threshold);
        out.field("surface_subdivisions",settings.surface_subdivisions);out.field("fit_radius",settings.fit_radius);out.field("max_bounces",settings.max_bounces);
        out.field("reflections_enabled",settings.reflections);
        out.field("particles",stats.particles);out.field("silicon_hits",stats.silicon_hits);
        out.field("mask_hits",stats.mask_hits);out.field("reflections",stats.reflections);
        out.field("escaped",stats.escaped);out.field("bounce_limit",stats.bounce_limit);
        out.field("interior_exits",transport.interior_exits);out.field("normal_mapped_hits",transport.normal_mapped_hits);
        out.field("mapping_failures",transport.mapping_failures);
        out.field("off_normal_mapped_hits",transport.off_normal_mapped_hits);out.field("max_mapping_offset_cells",transport.max_mapping_offset);
        out.field("pruned_si",stats.pruned_si);out.field("pruned_cl",stats.pruned_cl);
        out.field("periodic_orbits",transport.periodic_orbits);
        out.field("rank_deficient_fits",surface.rank_deficient_fits);out.field("undefined_normal_fits",surface.undefined_normal_fits);
        out.field("captures",reaction->captures);out.field("removed_si",reaction->removed);
        out.field("sicl4_products",reaction->products_sicl4);out.field("initial_si",initial_si);
        out.field("remaining_si",grid.silicon_inventory());out.field("surface_cl",grid.chlorine_inventory());
        out.field("relocated_cells",stats.relocated_cells);out.field("deleted_cells",stats.deleted_cells);
        out.field("fits",surface.fit_count);out.field("elapsed_seconds",elapsed_seconds);
    }
};
}
