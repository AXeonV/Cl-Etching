#include "planar/backend.hpp"
#include <functional>
#include <cfenv>
#include <algorithm>
#include <iostream>
#include <stdexcept>

using namespace etch;
using namespace etch::planar;

void check(bool condition, const std::string& message) {
  if (!condition)
    throw std::runtime_error(message);
}

void close(double actual, double expected, double tolerance = 1e-6) {
  check(std::abs(actual - expected) < tolerance,
    "expected " + std::to_string(expected) + ", got " + std::to_string(actual));
}

Geometry small_geometry() {
  return Geometry{40, 80, 20, 6, 16};
}

Settings small_settings() {
  Settings s{small_geometry()};
  s.max_bounces = 500;
  return s;
}

Source beam() {
  return {normal_direction(), constant(100)};
}

CellRemoval update_cell(CellUpdater& updater, const Grid& grid, int cell, Vec3 normal) {
  return updater.remove(
    cell, normal, {cell % grid.geometry.nx + 0.5, cell / grid.geometry.nx + 0.5});
}

void angular_functions() {
  close(cosine_yield(1), 1);
  close(cosine_yield(0), 0);
  close(enhanced_yield(1), 1);
  close(enhanced_yield(0), 0);
  close(cosine_yield(0.5), 0.5);
  close(enhanced_yield(0.5), 2);
}

// 每次反应后补回库存，使产额积分测试始终使用平面表面。
class FixedSurfaceYield : public YieldReaction {
public:
  using YieldReaction::YieldReaction;

  void apply(Sites& sites, double mu, double probability, Random& rng) override {
    YieldReaction::apply(sites, mu, probability, rng);
    sites[0] = 30;
  }
};

void yield_per_incident_ion() {
  const double pi = std::acos(-1.0);
  for (auto law : {cosine_yield, enhanced_yield}) {
    for (double degrees : {0.0, 30.0, 60.0, 80.0}) {
      for (bool fixed_capture : {false, true}) {
        Geometry g = small_geometry();
        g.opening = g.nx;
        double angle = degrees * pi / 180, mu = std::cos(angle);
        CaptureProbability probability =
          fixed_capture ? CaptureProbability([](double) { return 0.4; }) : cosine_capture;
        Simulation simulation(Settings{g},
          Source{
            [angle](Random&) { return Vec3{std::sin(angle), 0, std::cos(angle)}; }, constant(100)},
          std::make_unique<FixedSurfaceYield>(law), probability);
        simulation.advance(100000);
        close(simulation.reaction->removed / 100000.0, law(mu), 0.025);
        check(simulation.stats.silicon_hits == 100000, "one hit per fixed-plane incident ion");
        check(simulation.stats.deleted_cells == 0, "fixed-plane test does not change geometry");
        check(
          simulation.transport.diagnostics.interior_exits == 0, "fixed-plane rays are intercepted");
      }
    }
  }
  Grid grid(small_geometry());
  YieldReaction reaction(enhanced_yield);
  Random rng(8);
  int cell = grid.index(20, 20);
  grid.sites[cell][0] = 1;
  reaction.apply(grid.sites[cell], 0.5, 0.5, rng);
  check(grid.remaining(cell) == 0 && reaction.removed == 1, "removal is limited by cell inventory");
}

void source_distributions() {
  Random rng(42);
  double sigma = 10 * std::acos(-1.0) / 180 / 2.354820045;
  Source gaussian_source{gaussian_direction(sigma), constant(100)};
  Source cosine_source{cosine_direction(), constant(100)};
  double gaussian_x = 0, gaussian_y = 0, gaussian_x2 = 0, gaussian_y2 = 0;
  double cosine_z = 0, cosine_x2 = 0, cosine_y2 = 0;
  for (int i = 0; i < 100000; ++i) {
    Particle a = gaussian_source.emit(rng, small_geometry()),
             b = cosine_source.emit(rng, small_geometry());
    close(norm(a.direction), 1);
    close(norm(b.direction), 1);
    check(
      a.direction.z >= 0 && b.direction.z >= 0, "source directions lie in the downward hemisphere");
    check(a.position.x >= 0 && a.position.x < 40 && a.position.z == 0.25,
      "planar source positions remain two-dimensional");
    gaussian_x += a.direction.x;
    gaussian_y += a.direction.y;
    gaussian_x2 += a.direction.x * a.direction.x;
    gaussian_y2 += a.direction.y * a.direction.y;
    cosine_z += b.direction.z;
    cosine_x2 += b.direction.x * b.direction.x;
    cosine_y2 += b.direction.y * b.direction.y;
  }
  close(gaussian_x / 100000, 0, 0.001);
  close(gaussian_y / 100000, 0, 0.001);
  close(gaussian_x2 / 100000, sigma * sigma, 0.0002);
  close(gaussian_y2 / 100000, sigma * sigma, 0.0002);
  close(cosine_z / 100000, 2.0 / 3, 0.003);
  close(cosine_x2 / 100000, 0.25, 0.003);
  close(cosine_y2 / 100000, 0.25, 0.003);
}

void out_of_plane_collision() {
  Geometry g = small_geometry();
  g.opening = g.nx;
  Grid grid(g);
  Surface surface(grid);
  Transport transport(grid, surface);
  auto normal_hit = transport.next_hit({{20.3, 1}, {0, 0, 1}, 100});
  Particle particle{{20.3, 1}, {0, 0.8, 0.6}, 100};
  auto hit = transport.next_hit(particle);
  check(normal_hit && hit, "both projected rays reach the plane");
  close(norm(hit->position - normal_hit->position), 0);
  close(-dot(particle.direction, hit->normal), 0.6);
  Vec3 outgoing = reflected(particle.direction, hit->normal);
  close(outgoing.y, 0.8);
  close(outgoing.z, -0.6);
  close(norm(outgoing), 1);
  for (auto law : {cosine_yield, enhanced_yield}) {
    Simulation simulation(Settings{g},
      Source{[](Random&) { return Vec3{0, 0.8, 0.6}; }, constant(100)},
      std::make_unique<FixedSurfaceYield>(law));
    simulation.advance(100000);
    close(simulation.reaction->removed / 100000.0, law(0.6), 0.025);
    close(simulation.reaction->captures / 100000.0, 0.6, 0.01);
    check(simulation.stats.angles[53].primary == 100000,
      "collision angle uses all three direction components");
    check(simulation.stats.silicon_hits == 100000, "one collision per fixed-plane incident ray");
  }
}

void projected_horizontal_paths() {
  Geometry g = small_geometry();
  g.opening = g.nx;
  Grid grid(g);
  for (int z = 20; z < 50; ++z)
    for (int x = 12; x < 28; ++x) {
      int cell = grid.index(x, z);
      grid.material[cell] = Gas;
      grid.sites[cell] = {};
    }
  grid.refresh(20, 40, 80);
  Surface surface(grid);
  Transport transport(grid, surface);
  for (double dx : {0.1, 1e-12}) {
    Particle particle{{13, 30}, {dx, std::sqrt(1 - dx * dx), 0}, 100};
    auto wall = transport.next_hit(particle);
    check(wall.has_value(), "slow projected ray reaches the wall before its x period closes");
    close(wall->position.x, 28.5, 1e-5);
    close(wall->position.z, 30);
    Vec3 outgoing = reflected(particle.direction, wall->normal);
    close(outgoing.y, particle.direction.y);
    close(outgoing.x, -dx, 1e-15);
    auto opposite = transport.next_hit({wall->position, outgoing, 100});
    check(opposite.has_value(), "slow reflected projection reaches the opposite wall");
    close(opposite->position.x, 11.5, 1e-5);
  }
  for (Vec3 direction : {Vec3{0, 1, 0}, Vec3{0, -1, 0}, Vec3{0.1, std::sqrt(0.99), 0}}) {
    check(!transport.next_hit({{13, 1}, direction, 100}) &&
            transport.diagnostics.last_end == RayEnd::PeriodicOrbit,
      "unintercepted projected orbit terminates");
  }
}

void curve_geometry() {
  Curve flat{{0, 5}, {1, 0}, {0, -1}, 0, 0, 0};
  close(norm(flat.position(2) - Vec{2, 5}), 0);
  Curve parabola{{0, 5}, {1, 0}, {0, -1}, 0.2, 0, 0};
  close(norm(parabola.position(0.5) - Vec{0.5, 4.95}), 0);
  Vec n = parabola.normal(0.5);
  close(n.x, -0.2 / std::sqrt(1.04));
  close(n.z, -1 / std::sqrt(1.04));
}

void fitted_flat_surface() {
  Grid grid(small_geometry());
  Surface surface(grid);
  Vec normal = surface.curve(20, 20).normal(0);
  close(normal.x, 0);
  close(normal.z, -1);
  Transport transport(grid, surface);
  auto hit = transport.next_hit({{20.3, 1}, {0, 0, 1}, 100});
  check(hit.has_value(), "flat transport hit");
  close(hit->position.z, 20.5, 1e-5);
}

void vertical_wall_and_second_hit() {
  Grid grid(small_geometry());
  for (int z = 20; z < 50; ++z)
    for (int x = 12; x < 28; ++x) {
      int i = grid.index(x, z);
      grid.material[i] = Gas;
      grid.sites[i] = {};
    }
  grid.refresh(20, 40, 80);
  Surface surface(grid);
  Transport transport(grid, surface);
  Particle p{{25, 30}, unit(Vec3{1, 0, 2}), 100};
  auto wall = transport.next_hit(p);
  check(wall.has_value(), "wall collision");
  close(wall->position.x, 28.5, 1e-5);
  close(wall->normal.x, -1, 1e-5);
  p.direction = reflected(p.direction, wall->normal);
  p.position = wall->position;
  auto floor = transport.next_hit(p);
  check(floor.has_value(), "second collision after reflection");
  close(floor->position.z, 50.5, 0.2);
}

void v_groove_reflection_and_reverse_path() {
  Geometry g{100, 120, 20, 6, 100};
  Grid grid(g);
  for (int z = 20; z < 80; ++z)
    for (int x = 0; x < 100; ++x) {
      if (std::abs(x - 50) < 80 - z) {
        int cell = grid.index(x, z);
        grid.material[cell] = Gas;
        grid.sites[cell] = {};
      }
    }
  grid.refresh(50, 60, 120);
  Surface surface(grid);
  Transport transport(grid, surface);
  Particle ray{{60, 30}, {0, 0, 1}, 100};
  auto first = transport.next_hit(ray);
  check(first.has_value(), "V-groove first wall hit");
  close(first->position.x, 60);
  close(first->position.z, 71);
  ray.position = first->position;
  ray.direction = reflected(ray.direction, first->normal);
  close(ray.direction.x, -1);
  close(ray.direction.z, 0);
  auto second = transport.next_hit(ray);
  check(second.has_value(), "reflected ray hits opposite wall");
  close(second->position.x, 41);
  close(second->position.z, 71);
  Particle reverse{second->position, ray.direction * -1, 100};
  auto back = transport.next_hit(reverse);
  check(back.has_value(), "reversed ray returns to first wall");
  close(norm(back->position - first->position), 0);
  ray.position = second->position;
  ray.direction = reflected(ray.direction, second->normal);
  close(ray.direction.x, 0);
  close(ray.direction.z, -1);
  check(!transport.next_hit(ray) && transport.diagnostics.last_end == RayEnd::Escaped,
    "second reflection exits at top");
}

void reflection_geometry() {
  Vec3 d = unit(Vec3{1, 3, 2}), n = unit(Vec3{2, 0, -1}), r = reflected(d, n);
  close(norm(r), 1);
  close(dot(r, n), -dot(d, n));
  close(r.y, d.y);
  close(norm(reflected(r, n) - d), 0);
}

void surface_cell_mapping() {
  Grid grid(small_geometry());
  Surface surface(grid);
  Transport transport(grid, surface);
  Vec origin{20.5, 20.5};
  auto direct = transport.reaction_cell({20.7, 20.5}, {0, -1}, origin);
  check(direct && direct->x == 20 && direct->z == 20, "intersection cell takes precedence");
  auto gas = transport.reaction_cell({20.7, 19.7}, {0, -1}, origin);
  check(gas && gas->x == 20 && gas->z == 20, "gas intersection maps inward along normal");
  auto interior = transport.reaction_cell({20.7, 21.2}, {0, -1}, origin);
  check(interior && interior->x == 20 && interior->z == 20,
    "interior intersection maps outward along normal");
  auto closest = transport.reaction_cell({20.7, 19.7}, {1, 0}, {20.5, 19.7});
  check(closest && closest->z == 20 && transport.diagnostics.off_normal_mapped_hits == 1,
    "a normal that misses all cells maps to the closest surface cell to its line");
  auto unmapped = transport.reaction_cell({20.7, 10}, {1, 0}, {20.5, 10});
  check(!unmapped, "no arbitrary owner is substituted when the normal misses the surface");
}

void periodic_cell_mapping() {
  Geometry g = small_geometry();
  g.opening = g.nx;
  Grid grid(g);
  Surface surface(grid);
  Transport transport(grid, surface);
  auto cell = transport.reaction_cell({-0.2, 20.5}, {0, -1}, {-0.5, 20.5});
  check(cell && cell->x == g.nx - 1 && cell->z == 20, "periodic point maps to wrapped cell");
}

void chemistry_stoichiometry() {
  Geometry g = small_geometry();
  g.atoms_per_cell = 1;
  Grid grid(g);
  ChlorinationReaction reaction;
  Random rng(7);
  int cell = grid.index(20, 20);
  for (int state = 1; state <= 3; ++state) {
    reaction.apply(grid.sites[cell], 1, 1, rng);
    check(
      grid.sites[cell][state] == 1 && grid.remaining(cell) == 1, "single Cl advances one state");
  }
  reaction.apply(grid.sites[cell], 1, 1, rng);
  check(grid.remaining(cell) == 0 && reaction.removed == 1, "SiCl4 detachment");
  check(reaction.products_sicl4 == 1 && reaction.captures == 4, "four Cl per SiCl4");
}

void inventories_and_repeatability() {
  for (int model = 0; model < 2; ++model) {
    auto make = [model]() -> std::unique_ptr<Reaction> {
      if (model == 0)
        return std::make_unique<YieldReaction>();
      return std::make_unique<ChlorinationReaction>();
    };
    Simulation a(small_settings(), beam(), make()), b(small_settings(), beam(), make());
    a.advance(15000);
    b.advance(15000);
    check(a.grid.sites == b.grid.sites && a.grid.material == b.grid.material, "repeatable seed");
    check(a.grid.silicon_inventory() + a.reaction->removed + a.stats.detached_si == a.initial_si,
      "Si inventory");
    if (model == 1)
      check(a.reaction->captures ==
              a.grid.chlorine_inventory() + 4 * a.reaction->products_sicl4 + a.stats.detached_cl,
        "Cl inventory");
  }
}

void fit_cache_updates() {
  Grid grid(small_geometry());
  Surface surface(grid);
  surface.curve(20, 20);
  auto count = surface.fit_count;
  surface.curve(20, 20);
  check(surface.fit_count == count, "reuse unchanged fit");
  surface.candidates(20, 19);
  grid.sites[grid.index(20, 20)] = {};
  grid.remove_empty_cell(20, 20);
  surface.invalidate(20, 20);
  const auto& candidates = surface.candidates(20, 19);
  check(std::find(candidates.begin(), candidates.end(), grid.index(20, 20)) == candidates.end(),
    "remove stale broad-phase candidate");
  check(std::find(candidates.begin(), candidates.end(), grid.index(20, 21)) != candidates.end(),
    "index newly exposed surface cell");
  surface.curve(20, 21);
  check(surface.fit_count > count, "refit changed boundary");
}

void isolated_cell_has_no_centroid_normal() {
  Grid grid(small_geometry());
  int x = 20, z = 10, cell = grid.index(x, z);
  grid.material[cell] = Silicon;
  grid.sites[cell][0] = 12;
  grid.refresh(x, z, 2);
  Surface surface(grid);
  const Curve& curve = surface.curve(x, z);
  check(!curve.defined, "isolated cell does not receive an arbitrary normal");
  check(surface.undefined_normal_fits == 1, "record centroid degeneracy");
  check(grid.remaining(cell) == 12, "geometric diagnosis preserves isolated Si inventory");
}

void curve_normal_distance() {
  Curve flat{{0, 0}, {1, 0}, {0, 1}, 0, 0, 0};
  close(flat.project({2, 3}).distance, 3);
  close(flat.project({2, -3}).distance, -3);
  Curve parabola{{0, 0}, {1, 0}, {0, 1}, 1, 0, 0};
  auto projected = parabola.project({0, 2});
  close(std::abs(projected.u), std::sqrt(1.5));
  close(projected.distance, std::sqrt(1.75));
  close(dot(Vec{0, 2} - parabola.position(projected.u), unit(Vec{1, 2 * projected.u})), 0);
  close(parabola.project({1, 1}).distance, 0);
  close(parabola.project({0, 0.5}).u, 0);
  flat.a = 1e-18;
  close(flat.project({2, 3}).u, 2);
  close(flat.project({2, 3}).distance, 3);
}

void cell_removal_and_inventory_transfer() {
  for (int x : {20, 0}) {
    Geometry g = small_geometry();
    g.opening = g.nx;
    Grid grid(g);
    Surface surface(grid);
    CellUpdater updater(grid, surface);
    int original = grid.index(x, 20), protrusion = grid.index(x - 1, 19);
    grid.sites[original] = {1, 1, 1, 1};
    grid.material[protrusion] = Silicon;
    std::array<uint16_t, 4> inventory{10, 4, 5, 6};
    grid.sites[protrusion] = inventory;
    int intermediate = grid.index(x - 1, 20);
    Sites near_inventory{6, 7, 8, 9};
    grid.sites[intermediate] = near_inventory;
    grid.refresh(x, 20, 3);
    long long si = grid.silicon_inventory(), cl = grid.chlorine_inventory();
    surface.curve(x - 1, 20);
    auto fits = surface.fit_count;
    auto removal = update_cell(updater, grid, grid.index(x, 20), {0, 0, -1});
    check(removal.relocated && removal.cell == protrusion,
      "remove protruding neighbor rather than create a hole beneath it");
    check(grid.sites[original] == Sites{11, 5, 6, 7} &&
            grid.sites[intermediate] == near_inventory && grid.material[original] == Silicon,
      "move donor and preserve residual states without changing full intermediate cells");
    check(grid.material[protrusion] == Gas && grid.remaining(protrusion) == 0,
      "selected protrusion becomes gas");
    check(grid.silicon_inventory() == si && grid.chlorine_inventory() == cl,
      "all Si and Cl survive conservative relocation");
    surface.curve(x - 1, 20);
    check(surface.fit_count > fits, "invalidate fits around the actual removal location");
    for (int z = 18; z <= 22; ++z)
      for (int dx = -2; dx <= 2; ++dx) {
        int xx = x + dx, cell = grid.index(xx, z);
        bool boundary =
          grid.at(xx, z) != Gas && (grid.at(xx - 1, z) == Gas || grid.at(xx + 1, z) == Gas ||
                                     grid.at(xx, z - 1) == Gas || grid.at(xx, z + 1) == Gas);
        check(grid.boundary[cell] == boundary, "refresh the selected cell's neighborhood");
      }
  }
}

void flat_removal_and_fixed_mask() {
  Geometry g = small_geometry();
  g.opening = g.nx;
  Grid grid(g);
  Surface surface(grid);
  CellUpdater updater(grid, surface);
  grid.sites[grid.index(20, 20)] = {};
  auto removal = update_cell(updater, grid, grid.index(20, 20), {0, 0, -1});
  check(!removal.relocated && removal.cell == grid.index(20, 20),
    "flat surface advances at the underfilled cell");
  Grid masked(g);
  int mask = masked.index(19, 19);
  masked.material[mask] = Mask;
  masked.refresh(20, 20, 3);
  Surface masked_surface(masked);
  CellUpdater masked_updater(masked, masked_surface);
  masked.sites[masked.index(20, 20)] = {};
  update_cell(masked_updater, masked, masked.index(20, 20), {0, 0, -1});
  check(masked.material[mask] == Mask, "fixed mask is not a silicon removal candidate");
}

void connected_donor_selection() {
  Random selection_rng(17);
  for (int x : {8, 0}) {
    Grid grid({16, 24, 20, 0, 16});
    std::fill(grid.material.begin(), grid.material.end(), Gas);
    std::fill(grid.sites.begin(), grid.sites.end(), Sites{});
    int z = 8, bridge = grid.index(x, z), leaf = grid.index(x - 1, z);
    for (int dx : {-1, 0, 1})
      grid.material[grid.index(x + dx, z)] = Silicon;
    RemovalConnectivity<Grid> connectivity(grid);
    check(!connectivity.preserves_connection(bridge), "cutting a neck would detach a solid branch");
    std::vector<RemovalCandidate> candidates{{bridge, -1}, {leaf, 0}};
    check(connectivity.select(candidates, selection_rng).cell == leaf,
      "redirect removal to a connected leaf");
    check(connectivity.select(candidates, selection_rng).redirected &&
            connectivity.select(candidates, selection_rng).connected,
      "record topology redirection");
    candidates = {{bridge, 0}};
    check(connectivity.select(candidates, selection_rng).cell == bridge &&
            !connectivity.select(candidates, selection_rng).connected,
      "an unavoidable physical neck break still advances");
    for (int zz = 4; zz < z; ++zz)
      for (int dx : {-1, 1})
        grid.material[grid.index(x + dx, zz)] = Silicon;
    for (int dx = -1; dx <= 1; ++dx)
      grid.material[grid.index(x + dx, 4)] = Silicon;
    check(connectivity.preserves_connection(bridge),
      "distant alternate path avoids a false local disconnection, including x seam");
  }
  Grid grid({16, 24, 20, 0, 16});
  std::fill(grid.material.begin(), grid.material.end(), Gas);
  int original = grid.index(8, 8), attached = grid.index(7, 8), diagonal = grid.index(9, 7),
      mask = grid.index(9, 8);
  grid.material[original] = grid.material[attached] = grid.material[diagonal] = Silicon;
  grid.material[mask] = Mask;
  RemovalConnectivity<Grid> connectivity(grid);
  auto branch = connectivity.local_silicon_branch(original, {original, attached, diagonal, mask});
  check(branch.cells.size() == 2 && branch.contains(attached),
    "stock cannot cross a gas gap or fixed mask");
}

void collision_normal_at_symmetric_branch() {
  for (int sign : {-1, 1}) {
    Grid grid({20, 32, 24, 0, 20});
    std::fill(grid.material.begin(), grid.material.end(), Gas);
    std::fill(grid.sites.begin(), grid.sites.end(), Sites{});
    auto cell = [&](int x) { return grid.index(x, 8); };
    for (int x = 7; x <= 9; ++x) {
      grid.material[cell(x)] = Silicon;
      grid.sites[cell(x)] = {26, 0, 0, 0};
    }
    grid.sites[cell(8)] = {1, 1, 1, 1};
    grid.refresh(8, 8, 4);
    Surface surface(grid);
    CellUpdater updater(grid, surface);
    auto result = update_cell(updater, grid, cell(8), {double(sign), 0, 0});
    check(result.cell == cell(8 + sign),
      "collision normal orients removal when the local centroid is symmetric");
    check(updater.diagnostics.undefined_fit_removals == 1 &&
            updater.diagnostics.unavoidable_splits == 0,
      "singular local fit retains surface connection");
  }
}

void fragment_release_and_anchors() {
  Grid grid({20, 32, 24, 0, 20});
  std::fill(grid.material.begin(), grid.material.end(), Gas);
  std::fill(grid.sites.begin(), grid.sites.end(), Sites{});
  auto cell = [&](int z) { return grid.index(8, z); };
  for (int z = 5; z <= 11; ++z) {
    grid.material[cell(z)] = Silicon;
    grid.sites[cell(z)] = {6, 7, 8, 9};
  }
  grid.material[cell(12)] = Mask;
  grid.sites[cell(8)] = {};
  grid.refresh(8, 8, 8);
  Surface surface(grid);
  CellUpdater updater(grid, surface);
  for (int zz = 4; zz <= 9; ++zz)
    for (int xx = 7; xx <= 9; ++xx)
      surface.node(xx, zz);
  auto si = grid.silicon_inventory(), cl = grid.chlorine_inventory();
  auto removal = update_cell(updater, grid, cell(8), {0, 0, -1});
  check(removal.cell == cell(7) && removal.detached_cells == 2,
    "a real neck break releases the unanchored component");
  check(removal.detached_si == 60 && removal.detached_cl == 100,
    "released fragment states are accounted separately from reaction products");
  check(grid.silicon_inventory() + removal.detached_si == si, "Si closes after fracture");
  check(grid.chlorine_inventory() + removal.detached_cl == cl, "bound Cl closes after fracture");
  check(grid.material[cell(5)] == Gas && grid.material[cell(6)] == Gas &&
          grid.material[cell(12)] == Mask,
    "release only the unsupported branch");
  Surface fresh(grid);
  for (int zz = 4; zz <= 9; ++zz)
    for (int xx = 7; xx <= 9; ++xx) {
      auto cached = surface.node(xx, zz), expected = fresh.node(xx, zz);
      for (size_t derivative = 0; derivative < cached.size(); ++derivative)
        close(cached[derivative], expected[derivative], 1e-12);
    }
  RemovalConnectivity<Grid> connectivity(grid);
  check(
    connectivity.detached_components(cell(7)).empty(), "fragment release reaches a fixed point");
  grid.material[cell(5)] = Mask;
  for (int z : {6, 7})
    grid.material[cell(z)] = Silicon;
  grid.material[cell(8)] = Gas;
  check(connectivity.detached_components(cell(8)).empty(),
    "separate components attached to different fixed anchors both survive");
}

void conservative_update_defers_geometry() {
  Grid grid({20, 32, 12, 0, 20});
  int original = grid.index(10, 12);
  grid.sites[original] = {1, 1, 1, 1};
  Surface surface(grid);
  CellUpdater updater(grid, surface, 0.3, 17);
  surface.curve(10, 12);
  auto fits = surface.fit_count;
  auto before = grid.sites;
  auto material = grid.material;
  auto pending = update_cell(updater, grid, original, {0, 0, -1});
  check(pending.reaction_cells == 0 && pending.detached_cells == 0,
    "capacity shortage causes no geometric or inventory loss");
  check(grid.sites == before && grid.material == material && surface.fit_count == fits,
    "deferred updates preserve inventory and geometry caches");
  check(updater.diagnostics.deferred_updates == 1,
    "defer is recorded separately from a committed deletion");
  grid.sites[original] = {};
  auto committed = update_cell(updater, grid, original, {0, 0, -1});
  check(committed.reaction_cells == 1 && grid.material[original] == Gas,
    "empty cell still advances after a deferred update");
}

void unbiased_candidate_ties_and_capacity_support() {
  Grid grid({20, 32, 12, 0, 20});
  auto cell = [&](int x) { return grid.index(x, 12); };
  RemovalConnectivity<Grid> connectivity(grid);
  Random rng(23);
  int first = 0;
  for (int i = 0; i < 1000; ++i) {
    std::vector<RemovalCandidate> candidates{{cell(9), 0, 1}, {cell(9), 0, 1}, {cell(11), 0, 1}};
    first += connectivity.select(candidates, rng).cell == cell(9);
  }
  check(first > 400 && first < 600,
    "periodic aliases do not bias otherwise equivalent deletion choices");
  std::vector<RemovalCandidate> candidates{{cell(9), 0, 4}, {cell(11), 0, 1}};
  auto before_rng = rng;
  check(connectivity.select(candidates, rng).cell == cell(11) && rng == before_rng,
    "minimize transport among equal geometric scores without unnecessary RNG draws");
  grid.sites[cell(10)] = {1, 1, 1, 1};
  grid.sites[cell(12)] = {25, 0, 0, 0};
  Surface surface(grid);
  CellUpdater local(grid, surface, 0.3, 17);
  auto pending = update_cell(local, grid, cell(10), {0, 0, -1});
  check(pending.reaction_cells == 0, "strictly local remapping defers when its own volume is full");
}

void connectivity_matches_global_components() {
  Grid grid({7, 9, 5, 0, 7});
  Random rng(77);
  RemovalConnectivity<Grid> connectivity(grid);
  auto components = [&](int removed) {
    std::vector<bool> seen(grid.material.size());
    int count = 0;
    for (size_t seed = 0; seed < grid.material.size(); ++seed)
      if (int(seed) != removed && grid.material[seed] != Gas && !seen[seed]) {
        ++count;
        std::vector<int> queue{int(seed)};
        seen[seed] = true;
        for (size_t i = 0; i < queue.size(); ++i)
          for (int neighbor : grid.face_neighbors(queue[i])) {
            if (neighbor >= 0 && neighbor != removed && grid.material[neighbor] != Gas &&
                !seen[neighbor]) {
              seen[neighbor] = true;
              queue.push_back(neighbor);
            }
          }
      }
    return count;
  };
  for (int sample = 0; sample < 500; ++sample) {
    for (auto& material : grid.material) {
      double value = uniform(rng);
      material = value < 0.4 ? Silicon : value < 0.45 ? Mask : Gas;
    }
    int before = components(-1), checked = 0;
    for (size_t n = 0; n < grid.material.size() && checked < 5; ++n) {
      int cell = int((n + sample) % grid.material.size());
      if (grid.material[cell] != Silicon)
        continue;
      check(connectivity.preserves_connection(cell) == (components(cell) <= before),
        "multi-source decision equals full component enumeration");
      ++checked;
    }
  }
}

void collision_projection_and_periodic_image() {
  std::array<int, 2> chosen{};
  for (int side = 0; side < 2; ++side) {
    Grid grid({40, 80, 10, 0, 40});
    for (int z = 0; z < 80; ++z)
      for (int x = 0; x < 40; ++x) {
        int cell = grid.index(x, z), top = std::max(0, 32 - (x - 20) * (x - 20));
        grid.material[cell] = z >= top ? Silicon : Gas;
        grid.sites[cell] = z >= top ? Sites{30, 0, 0, 0} : Sites{};
      }
    for (int x : {19, 21}) {
      int cell = grid.index(x, 30);
      grid.material[cell] = Silicon;
      grid.sites[cell] = {30, 0, 0, 0};
    }
    int original = grid.index(20, 32);
    grid.sites[original] = {};
    grid.refresh(20, 40, 80);
    Surface surface(grid);
    Curve curve = surface.curve(20, 32);
    Vec point = curve.position(side == 0 ? -0.3 : 0.3);
    CellUpdater updater(grid, surface);
    chosen[side] = updater.remove(original, {0, 0, -1}, point).cell;
  }
  check(chosen[0] != chosen[1] && chosen[0] % 40 + chosen[1] % 40 == 40,
    "the actual contact projection selects reflected counterparts on a curved floor");
  Geometry g{40, 80, 20, 0, 40};
  Grid a(g), b(g);
  for (auto* grid : {&a, &b}) {
    int tip = grid->index(-1, 19);
    grid->material[tip] = Silicon;
    grid->sites[tip] = {30, 0, 0, 0};
    grid->sites[grid->index(0, 20)] = {};
    grid->refresh(0, 20, 4);
  }
  Surface sa(a), sb(b);
  CellUpdater ua(a, sa), ub(b, sb);
  ua.remove(a.index(0, 20), {0, 0, -1}, {0.25, 20.5});
  ub.remove(b.index(0, 20), {0, 0, -1}, {40.25, 20.5});
  check(a.material == b.material && a.sites == b.sites,
    "equivalent periodic contacts have identical updates");
}

void undefined_normal_removal() {
  for (bool symmetric_neighbors : {false, true}) {
    Grid grid(small_geometry());
    int x = 20, z = 10, original = grid.index(x, z);
    grid.material[original] = Silicon;
    grid.sites[original] = {1, 1, 1, 1};
    if (symmetric_neighbors)
      for (int dx : {-1, 1}) {
        int cell = grid.index(x + dx, z);
        grid.material[cell] = Silicon;
        grid.sites[cell] = {8, 7, 6, 5};
      }
    grid.refresh(x, z, 3);
    Surface surface(grid);
    CellUpdater updater(grid, surface);
    check(!surface.curve(x, z).defined, "symmetric neighborhood has no centroid normal");
    auto si = grid.silicon_inventory(), cl = grid.chlorine_inventory();
    std::feclearexcept(FE_ALL_EXCEPT);
    auto result = update_cell(updater, grid, original, {0, 0, -1});
    check(!std::fetestexcept(FE_INVALID | FE_DIVBYZERO),
      "undefined fit is not projected or normalized");
    if (symmetric_neighbors) {
      check(result.cell != original && result.relocated,
        "a symmetric connected structure uses the collision normal and topology");
      check(updater.diagnostics.unavoidable_splits == 0,
        "undefined centroid does not force a neck break");
    } else {
      check(
        result.cell == original && grid.material[original] == Gas && grid.remaining(original) == 0,
        "a true singleton is cleared");
      check(result.reaction_cells == 0 && result.detached_si == 4 && result.detached_cl == 6,
        "a singleton leaves as solid fragment inventory rather than chemical residual");
    }
    check(grid.silicon_inventory() + result.detached_si == si &&
            grid.chlorine_inventory() + result.detached_cl == cl,
      "Si and bound Cl close for connected and isolated updates");
    if (symmetric_neighbors)
      check(grid.remaining(original) >= 26 && grid.remaining(original) <= 30,
        "the surviving branch retains donor inventory plus conserved residual states");
  }
}

int main() {
  std::cout << std::unitbuf;
  std::vector<std::pair<std::string, std::function<void()>>> tests = {
    {"actual collision projection and periodic image", collision_projection_and_periodic_image},
    {"multi-source connectivity equivalence", connectivity_matches_global_components},
    {"unbiased choices and capacity support", unbiased_candidate_ties_and_capacity_support},
    {"conservative deferred geometry", conservative_update_defers_geometry},
    {"fragment inventory and fixed anchors", fragment_release_and_anchors},
    {"collision normal at symmetric branch", collision_normal_at_symmetric_branch},
    {"connected inventory donor selection", connected_donor_selection},
    {"undefined normal removal", undefined_normal_removal},
    {"centroid degeneracy", isolated_cell_has_no_centroid_normal},
    {"normal distance to quadratic curve", curve_normal_distance},
    {"cell removal and inventory transfer", cell_removal_and_inventory_transfer},
    {"flat removal and fixed mask", flat_removal_and_fixed_mask},
    {"angular functions", angular_functions},
    {"yield per incident ion including capture and reflection", yield_per_incident_ion},
    {"source distribution moments", source_distributions},
    {"three-component collision and capture", out_of_plane_collision},
    {"projected horizontal periods", projected_horizontal_paths},
    {"local quadratic geometry", curve_geometry}, {"fitted flat surface", fitted_flat_surface},
    {"vertical fit and secondary collision", vertical_wall_and_second_hit},
    {"V-groove reflection and reverse path", v_groove_reflection_and_reverse_path},
    {"specular reflection", reflection_geometry},
    {"single Cl stoichiometry", chemistry_stoichiometry},
    {"inventory and seed", inventories_and_repeatability}, {"fit cache update", fit_cache_updates},
    {"reaction-cell normal mapping", surface_cell_mapping},
    {"periodic cell mapping", periodic_cell_mapping}};
  int failed = 0;
  for (auto& [name, test] : tests) {
    try {
      test();
      std::cout << "PASS " << name << '\n';
    } catch (const std::exception& e) {
      ++failed;
      std::cerr << "FAIL " << name << ": " << e.what() << '\n';
    }
  }
  std::cout << tests.size() - failed << "/" << tests.size() << " test groups passed\n";
  return failed ? 1 : 0;
}
