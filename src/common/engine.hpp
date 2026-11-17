#pragma once
#include "reaction.hpp"
#include "events.hpp"
#include "output.hpp"
#include <filesystem>
#include <memory>
#include <algorithm>

namespace etch {
// 逐个追踪粒子，反应、反射和元胞更新共享同一碰撞结果。
template <class Backend> class Engine {
public:
  using Settings = RunSettings<typename Backend::Geometry>;
  Settings settings;
  typename Backend::Grid grid;
  typename Backend::Surface surface;
  typename Backend::Transport transport;
  typename Backend::CellUpdater cell_updater;
  typename Backend::Source source;
  CaptureProbability capture;
  std::unique_ptr<Reaction> reaction;
  // 发射、反应和元胞更新分别使用独立随机流，避免额外抽样改变入射序列。
  Random source_rng, event_rng;
  Statistics stats;
  std::vector<uint32_t> primary_hits, secondary_hits;
  long long initial_si;

  Engine(Settings s, typename Backend::Source src, std::unique_ptr<Reaction> rx,
    CaptureProbability p = cosine_capture)
      : settings(s), grid(s.geometry), surface(grid, s.fit_radius, s.surface_subdivisions),
        transport(grid, surface), cell_updater(grid, surface, s.smoothing_weight, s.seed + 2000003),
        source(std::move(src)), capture(std::move(p)), reaction(std::move(rx)), source_rng(s.seed),
        event_rng(s.seed + 1000003), primary_hits(grid.material.size()),
        secondary_hits(grid.material.size()), initial_si(grid.silicon_inventory()) {}

  void advance(long long particle_count) {
    for (long long particle_index = 0; particle_index < particle_count; ++particle_index) {
      auto particle = source.emit(source_rng, grid.geometry);
      ++stats.particles;

      for (int bounce = 0;; ++bounce) {
        auto hit = transport.next_hit(particle);
        if (!hit) {
          if (transport.diagnostics.last_end == RayEnd::Escaped)
            ++stats.escaped;
          break;
        }
        int cell = hit->cell;
        // 2D3V 仍用完整三分量方向计算入射角，法线指向气体侧。
        double mu = -dot(particle.direction, hit->normal);
        if (grid.material[cell] == Mask)
          ++stats.mask_hits;
        else {
          ++stats.silicon_hits;
          int angle_bin =
            std::min(89, int(std::acos(std::clamp(mu, 0.0, 1.0)) * 180 / std::acos(-1.0)));
          auto& bin = stats.angles[angle_bin];
          ++(bounce == 0 ? bin.primary : bin.secondary);
          (bounce == 0 ? primary_hits : secondary_hits)[cell]++;

          double capture_probability = capture(mu);
          if (uniform(event_rng) < capture_probability) {
            long long removed_before = reaction->removed;
            reaction->apply(grid.sites[cell], mu, capture_probability, event_rng);
            ++bin.captures;
            bin.removed += reaction->removed - removed_before;

            // 反应先改变库存；低于阈值后才尝试几何更新或延后删格。
            if (grid.remaining(cell) < grid.geometry.deletion_threshold) {
              auto removal = cell_updater.remove(cell, hit->normal, hit->position);
              stats.relocated_cells += removal.relocated;
              stats.deleted_cells += removal.reaction_cells;
              stats.detached_cells += removal.detached_cells;
              stats.detached_si += removal.detached_si;
              stats.detached_cl += removal.detached_cl;
              stats.redistributed_si += removal.redistributed_si;
              stats.redistributed_cl += removal.redistributed_cl;
            }
            // 捕获后结束该粒子的行程，即使本次没有删除材料格。
            break;
          }
        }

        // 固定掩膜或未捕获的粒子继续镜面反射，能量保持原值。
        if (!settings.reflections)
          break;
        if (settings.max_bounces > 0 && bounce == settings.max_bounces) {
          ++stats.bounce_limit;
          break;
        }
        particle.position = hit->position;
        grid.wrap_position(particle.position);
        particle.direction = reflected(particle.direction, hit->normal);
        ++stats.reflections;
      }
    }
  }

  void save(const std::string& directory, double elapsed_seconds) {
    std::filesystem::create_directories(directory);
    write_binary(directory + "/material.bin", grid.material);
    write_binary(directory + "/sites.bin", grid.sites);
    write_binary(directory + "/primary_hits.bin", primary_hits);
    write_binary(directory + "/secondary_hits.bin", secondary_hits);

    Backend::write_surface(surface, grid, directory + "/surface.csv");
    Backend::write_field(surface, grid, directory + "/surface_nodes.bin");

    std::ofstream angles(directory + "/collision_angles.csv");
    angles << "angle_lower_deg,angle_upper_deg,primary_hits,secondary_hits,captures,removed_si\n";
    for (size_t i = 0; i < stats.angles.size(); ++i) {
      const auto& bin = stats.angles[i];
      angles << i << ',' << i + 1 << ',' << bin.primary << ',' << bin.secondary << ','
             << bin.captures << ',' << bin.removed << '\n';
    }

    Summary out(directory + "/summary.json");
    out.field("model", reaction->name());
    out.field("surface_representation", std::string("hermite"));
    out.field("direction_components", 3);
    Backend::write_geometry(out, grid.geometry);
    out.field("atoms_per_cell", grid.geometry.atoms_per_cell);
    out.field("seed", settings.seed);
    out.field("deletion_threshold", grid.geometry.deletion_threshold);
    out.field("surface_subdivisions", settings.surface_subdivisions);
    out.field("fit_radius", settings.fit_radius);
    out.field("max_bounces", settings.max_bounces);
    out.field("reflections_enabled", settings.reflections);

    out.field("particles", stats.particles);
    out.field("silicon_hits", stats.silicon_hits);
    out.field("mask_hits", stats.mask_hits);
    out.field("reflections", stats.reflections);
    out.field("escaped", stats.escaped);
    out.field("bounce_limit", stats.bounce_limit);
    out.field("interior_exits", transport.diagnostics.interior_exits);
    out.field("normal_mapped_hits", transport.diagnostics.normal_mapped_hits);
    out.field("mapping_failures", transport.diagnostics.mapping_failures);
    out.field("off_normal_mapped_hits", transport.diagnostics.off_normal_mapped_hits);
    out.field("max_mapping_offset_cells", transport.diagnostics.max_mapping_offset);

    out.field("detached_si", stats.detached_si);
    out.field("detached_cl", stats.detached_cl);
    out.field("detached_cells", stats.detached_cells);
    out.field("detached_material", std::string("released_from_static_substrate"));
    out.field("inventory_update", std::string("conservative"));
    out.field("update_seed", settings.seed + 2000003);
    out.field("smoothing_weight", settings.smoothing_weight);
    out.field("inventory_remap", std::string("slots"));
    out.field("update_geometry", std::string("oriented_quadratic"));
    out.field("update_normal", std::string("projected-fit"));
    out.field("reoriented_fits", cell_updater.diagnostics.reoriented_fits);
    out.field("tangent_fallbacks", cell_updater.diagnostics.tangent_fallbacks);
    out.field("redistributed_si", stats.redistributed_si);
    out.field("redistributed_cl", stats.redistributed_cl);
    out.field("deferred_updates", cell_updater.diagnostics.deferred_updates);
    out.field("periodic_orbits", transport.diagnostics.periodic_orbits);
    out.field("rank_deficient_fits", surface.rank_deficient_fits);
    out.field("undefined_normal_fits", surface.undefined_normal_fits);

    // 化学去除与固体碎片分别记账，库存核对使用两者的总和。
    out.field("captures", reaction->captures);
    out.field("removed_si", reaction->removed);
    out.field("sicl4_products", reaction->products_sicl4);
    out.field("initial_si", initial_si);
    out.field("remaining_si", grid.silicon_inventory());
    out.field("surface_cl", grid.chlorine_inventory());
    out.field("relocated_cells", stats.relocated_cells);
    out.field("deleted_cells", stats.deleted_cells);
    out.field("cell_update_method", std::string("local_connected_transfer"));
    out.field("undefined_fit_removals", cell_updater.diagnostics.undefined_fit_removals);
    out.field("disconnected_donors", cell_updater.diagnostics.disconnected_donors);
    out.field("redirected_removals", cell_updater.diagnostics.redirected_removals);
    out.field("unavoidable_splits", cell_updater.diagnostics.unavoidable_splits);
    out.field("fits", surface.fit_count);
    out.field("elapsed_seconds", elapsed_seconds);
  }
};
}
