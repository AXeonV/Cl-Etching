#include "simulation.hpp"
#include <functional>
#include <algorithm>
#include <iostream>
#include <stdexcept>

using namespace etch;
void check(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}
void close(double actual, double expected, double tolerance = 1e-6) {
    check(std::abs(actual - expected) < tolerance, "expected " + std::to_string(expected)
          + ", got " + std::to_string(actual));
}
Geometry small_geometry() {
    return Geometry{40, 80, 20, 6, 16};
}
Settings small_settings() { Settings s{small_geometry()}; s.max_bounces = 500; return s; }
Source beam() { return {constant(0), constant(100)}; }

void angular_functions() {
    close(cosine_yield(1), 1); close(cosine_yield(0), 0);
    close(enhanced_yield(1), 1); close(enhanced_yield(0), 0);
    close(cosine_yield(0.5), 0.5); close(enhanced_yield(0.5), 2);
}
// Restore inventory after each event so this integration test keeps a flat surface.
class FixedSurfaceYield : public YieldReaction {
public:
    using YieldReaction::YieldReaction;
    void apply(Grid& grid, int cell, double mu, double probability, Random& rng) override {
        YieldReaction::apply(grid, cell, mu, probability, rng);
        grid.sites[cell][0] = grid.geometry.atoms_per_cell;
    }
};
void yield_per_incident_ion() {
    const double pi = std::acos(-1.0);
    for (auto law : {cosine_yield, enhanced_yield}) {
        for (double degrees : {0.0, 30.0, 60.0, 80.0}) {
            for (bool fixed_capture : {false, true}) {
                Geometry g = small_geometry(); g.opening = g.nx;
                double angle = degrees*pi/180, mu = std::cos(angle);
                CaptureProbability probability = fixed_capture
                    ? CaptureProbability([](double) { return 0.4; }) : cosine_capture;
                Simulation simulation(Settings{g}, Source{constant(angle),constant(100)},
                                      std::make_unique<FixedSurfaceYield>(law), probability);
                simulation.advance(100000);
                close(simulation.reaction->removed / 100000.0, law(mu), 0.025);
                check(simulation.stats.silicon_hits == 100000, "one hit per fixed-plane incident ion");
                check(simulation.stats.deleted_cells == 0, "fixed-plane test does not change geometry");
                check(simulation.transport.interior_exits == 0, "fixed-plane rays are intercepted");
            }
        }
    }
    Grid grid(small_geometry()); YieldReaction reaction(enhanced_yield); Random rng(8);
    int cell = grid.index(20,20);
    grid.sites[cell][0] = 1;
    reaction.apply(grid, cell, 0.5, 0.5, rng);
    check(grid.remaining(cell) == 0 && reaction.removed == 1, "removal is limited by cell inventory");
}
void source_distributions() {
    Random rng(42);
    double sigma = 10 * std::acos(-1.0) / 180 / 2.354820045;
    auto angle = gaussian(0, sigma, -1.5, 1.5);
    auto cosine = cosine_angle();
    double sum = 0, sum2 = 0, cosine_sin2 = 0;
    for (int i = 0; i < 100000; ++i) {
        double a = angle(rng); sum += a; sum2 += a * a;
        double b = cosine(rng); cosine_sin2 += std::pow(std::sin(b), 2);
    }
    close(sum / 100000, 0, 0.001);
    close(sum2 / 100000, sigma * sigma, 0.0001);
    close(cosine_sin2 / 100000, 1.0 / 3, 0.005);
}
void curve_intersection() {
    Curve flat{{0, 5}, {1, 0}, {0, -1}, 0, 0, 0};
    auto hit = flat.intersection({0, 0}, {0, 1}, 0, 10);
    check(hit.has_value(), "flat curve hit"); close(*hit, 5);
    check(!flat.intersection({0, 0}, {0, -1}, 0, 10), "backward root rejected");
    Curve parabola{{0, 5}, {1, 0}, {0, -1}, 0.2, 0, 0};
    auto curved_hit = parabola.intersection({0.5, 0}, {0, 1}, 0, 10);
    check(curved_hit.has_value(), "quadratic hit"); close(*curved_hit, 4.95);
    Vec n = parabola.normal(0.5);
    close(n.x, -0.2 / std::sqrt(1.04)); close(n.z, -1 / std::sqrt(1.04));
}
void fitted_flat_surface() {
    Grid grid(small_geometry()); Surface surface(grid);
    Vec normal = surface.curve(20, 20).normal(0);
    close(normal.x, 0); close(normal.z, -1);
    Transport transport(grid, surface);
    auto hit = transport.next_hit({{20.3, 1}, {0, 1}, 100});
    check(hit.has_value(), "flat transport hit"); close(hit->position.z, 20.5, 1e-5);
}
void vertical_wall_and_second_hit() {
    Grid grid(small_geometry());
    for (int z = 20; z < 50; ++z) for (int x = 12; x < 28; ++x) {
        int i = grid.index(x, z); grid.material[i] = Gas; grid.sites[i] = {};
    }
    grid.refresh(20, 40, 80);
    Surface surface(grid); Transport transport(grid, surface);
    Particle p{{25, 30}, unit({1, 2}), 100};
    auto wall = transport.next_hit(p);
    check(wall.has_value(), "wall collision");
    close(wall->position.x, 28.5, 1e-5); close(wall->normal.x, -1, 1e-5);
    p.direction = reflected(p.direction, wall->normal);
    p.position = wall->position;
    auto floor = transport.next_hit(p);
    check(floor.has_value(), "second collision after reflection");
    close(floor->position.z, 50.5, 0.2);
}
void v_groove_reflection_and_reverse_path() {
    Geometry g{100,120,20,6,100}; Grid grid(g);
    for (int z=20;z<80;++z) for(int x=0;x<100;++x) {
        if (std::abs(x-50)<80-z) {
            int cell=grid.index(x,z);grid.material[cell]=Gas;grid.sites[cell]={};
        }
    }
    grid.refresh(50,60,120);
    Surface surface(grid);Transport transport(grid,surface);
    Particle ray{{60,30},{0,1},100};
    auto first=transport.next_hit(ray);check(first.has_value(),"V-groove first wall hit");
    close(first->position.x,60);close(first->position.z,71);
    ray.position=first->position;ray.direction=reflected(ray.direction,first->normal);
    close(ray.direction.x,-1);close(ray.direction.z,0);
    auto second=transport.next_hit(ray);check(second.has_value(),"reflected ray hits opposite wall");
    close(second->position.x,41);close(second->position.z,71);
    Particle reverse{second->position,ray.direction*-1,100};
    auto back=transport.next_hit(reverse);check(back.has_value(),"reversed ray returns to first wall");
    close(norm(back->position-first->position),0);
    ray.position=second->position;ray.direction=reflected(ray.direction,second->normal);
    close(ray.direction.x,0);close(ray.direction.z,-1);
    check(!transport.next_hit(ray)&&transport.last_end==RayEnd::Escaped,"second reflection exits at top");
}

void reflection_geometry() {
    Vec d = unit({1, 2}), n{0, -1}, r = reflected(d, n);
    close(norm(r), 1); close(dot(r, n), -dot(d, n)); close(r.x, d.x);
}
void surface_cell_mapping() {
    Grid grid(small_geometry()); Surface surface(grid); Transport transport(grid, surface);
    Vec origin{20.5, 20.5};
    auto direct = transport.reaction_cell({20.7, 20.5}, {0, -1}, origin);
    check(direct && direct->x == 20 && direct->z == 20, "intersection cell takes precedence");
    auto gas = transport.reaction_cell({20.7, 19.7}, {0, -1}, origin);
    check(gas && gas->x == 20 && gas->z == 20, "gas intersection maps inward along normal");
    auto interior = transport.reaction_cell({20.7, 21.2}, {0, -1}, origin);
    check(interior && interior->x == 20 && interior->z == 20, "interior intersection maps outward along normal");
    auto unmapped = transport.reaction_cell({20.7, 10}, {1, 0}, {20.5, 10});
    check(!unmapped, "no arbitrary owner is substituted when the normal misses the surface");
}
void periodic_cell_mapping() {
    Geometry g = small_geometry(); g.opening = g.nx;
    Grid grid(g); Surface surface(grid); Transport transport(grid, surface);
    auto cell = transport.reaction_cell({-0.2, 20.5}, {0, -1}, {-0.5, 20.5});
    check(cell && cell->x == g.nx - 1 && cell->z == 20, "periodic point maps to wrapped cell");
}
void curve_support_is_set_by_the_tracer() {
    Curve flat{{0, 5}, {1, 0}, {0, -1}, 0, 0, 0};
    auto hit = flat.intersection({2, 0}, {0, 1}, 0, 10);
    check(hit.has_value(), "valid roots are not discarded by an extra one-cell clipping box");
    close(*hit, 5);
}
void chemistry_stoichiometry() {
    Geometry g = small_geometry(); g.atoms_per_cell = 1;
    Grid grid(g); ChlorinationReaction reaction; Random rng(7);
    int cell = grid.index(20, 20);
    for (int state = 1; state <= 3; ++state) {
        reaction.apply(grid, cell, 1, 1, rng);
        check(grid.sites[cell][state] == 1 && grid.remaining(cell) == 1, "single Cl advances one state");
    }
    reaction.apply(grid, cell, 1, 1, rng);
    check(grid.remaining(cell) == 0 && reaction.removed == 1, "SiCl4 detachment");
    check(reaction.products_sicl4 == 1 && reaction.captures == 4, "four Cl per SiCl4");
}
void inventories_and_repeatability() {
    for (int model = 0; model < 2; ++model) {
        auto make = [model]() -> std::unique_ptr<Reaction> {
            if (model == 0) return std::make_unique<YieldReaction>();
            return std::make_unique<ChlorinationReaction>();
        };
        Simulation a(small_settings(), beam(), make()), b(small_settings(), beam(), make());
        a.advance(15000); b.advance(15000);
        check(a.grid.sites == b.grid.sites && a.grid.material == b.grid.material, "repeatable seed");
        check(a.grid.silicon_inventory() + a.reaction->removed + a.stats.pruned_si == a.initial_si, "Si inventory");
        if (model == 1) check(a.reaction->captures == a.grid.chlorine_inventory()
            + 4 * a.reaction->products_sicl4 + a.stats.pruned_cl, "Cl inventory");
    }
}
void fit_cache_updates() {
    Grid grid(small_geometry()); Surface surface(grid);
    surface.curve(20, 20); auto count = surface.fit_count;
    surface.curve(20, 20); check(surface.fit_count == count, "reuse unchanged fit");
    surface.candidates(20,19);
    grid.sites[grid.index(20, 20)] = {};
    grid.remove_empty_cell(20, 20); surface.invalidate(20, 20);
    const auto& candidates=surface.candidates(20,19);
    check(std::find(candidates.begin(),candidates.end(),grid.index(20,20))==candidates.end(), "remove stale broad-phase candidate");
    check(std::find(candidates.begin(),candidates.end(),grid.index(20,21))!=candidates.end(), "index newly exposed surface cell");
    surface.curve(20, 21); check(surface.fit_count > count, "refit changed boundary");
}
void isolated_cell_has_no_centroid_normal() {
    Grid grid(small_geometry());
    int x=20,z=10,cell=grid.index(x,z);
    grid.material[cell]=Silicon;grid.sites[cell][0]=12;grid.refresh(x,z,2);
    Surface surface(grid);
    const Curve& curve=surface.curve(x,z);
    check(!curve.defined,"isolated cell does not receive an arbitrary normal");
    check(surface.undefined_normal_fits==1,"record centroid degeneracy");
    check(grid.remaining(cell)==12,"geometric diagnosis preserves isolated Si inventory");
}

void curve_normal_distance() {
    Curve flat{{0, 0}, {1, 0}, {0, 1}, 0, 0, 0};
    close(flat.project({2, 3}).distance, 3);
    close(flat.project({2, -3}).distance, -3);
    Curve parabola{{0, 0}, {1, 0}, {0, 1}, 1, 0, 0};
    auto projected = parabola.project({0, 2});
    close(std::abs(projected.u), std::sqrt(1.5));
    close(projected.distance, std::sqrt(1.75));
    close(dot(Vec{0, 2}-parabola.position(projected.u),
              unit(Vec{1, 2*projected.u})), 0);
    close(parabola.project({1, 1}).distance, 0);
    close(parabola.project({0, 0.5}).u, 0);
    flat.a = 1e-18;
    close(flat.project({2, 3}).u, 2);
    close(flat.project({2, 3}).distance, 3);
}
void cell_removal_and_inventory_transfer() {
    for (int x : {20, 0}) {
        Geometry g = small_geometry(); g.opening = g.nx;
        Grid grid(g); Surface surface(grid); CellUpdater updater(grid, surface);
        int original = grid.index(x, 20), protrusion = grid.index(x-1, 19);
        grid.sites[original] = {1, 1, 1, 1};
        grid.material[protrusion] = Silicon;
        std::array<uint16_t, 4> inventory{10, 4, 5, 6};
        grid.sites[protrusion] = inventory;
        grid.refresh(x, 20, 3);
        long long si = grid.silicon_inventory(), cl = grid.chlorine_inventory();
        surface.curve(x-1, 20); auto fits = surface.fit_count;
        auto removal = updater.remove(x, 20);
        check(removal.relocated && grid.index(removal.x, removal.z) == protrusion,
              "remove protruding neighbor rather than create a hole beneath it");
        check(grid.sites[original] == inventory && grid.material[original] == Silicon,
              "transfer all four chlorination states into the underfilled cell");
        check(grid.material[protrusion] == Gas && grid.remaining(protrusion) == 0,
              "selected protrusion becomes gas");
        check(grid.silicon_inventory()+4 == si && grid.chlorine_inventory()+6 == cl,
              "only the original subthreshold residual leaves the inventory");
        surface.curve(x-1, 20);
        check(surface.fit_count > fits, "invalidate fits around the actual removal location");
        for (int z = 18; z <= 22; ++z) for (int dx = -2; dx <= 2; ++dx) {
            int xx = x+dx, cell = grid.index(xx,z);
            bool boundary = grid.at(xx,z) != Gas && (grid.at(xx-1,z)==Gas || grid.at(xx+1,z)==Gas
                || grid.at(xx,z-1)==Gas || grid.at(xx,z+1)==Gas);
            check(grid.boundary[cell] == boundary, "refresh the selected cell's neighborhood");
        }
    }
}
void flat_removal_and_fixed_mask() {
    Geometry g = small_geometry(); g.opening = g.nx;
    Grid grid(g); Surface surface(grid); CellUpdater updater(grid, surface);
    grid.sites[grid.index(20,20)] = {};
    auto removal = updater.remove(20,20);
    check(!removal.relocated && removal.x == 20 && removal.z == 20,
          "flat surface advances at the underfilled cell");
    Grid masked(g); int mask = masked.index(19,19);
    masked.material[mask] = Mask; masked.refresh(20,20,3);
    Surface masked_surface(masked); CellUpdater masked_updater(masked,masked_surface);
    masked.sites[masked.index(20,20)] = {};
    masked_updater.remove(20,20);
    check(masked.material[mask] == Mask, "fixed mask is not a silicon removal candidate");
}

int main() {
    std::cout << std::unitbuf;
    std::vector<std::pair<std::string, std::function<void()>>> tests = {
        {"centroid degeneracy", isolated_cell_has_no_centroid_normal},
        {"normal distance to quadratic curve", curve_normal_distance},
        {"cell removal and inventory transfer", cell_removal_and_inventory_transfer},
        {"flat removal and fixed mask", flat_removal_and_fixed_mask},
        {"angular functions", angular_functions}, {"yield per incident ion including capture and reflection", yield_per_incident_ion},
        {"source distribution moments", source_distributions}, {"curve intersection and derivative", curve_intersection},
        {"fitted flat surface", fitted_flat_surface}, {"vertical fit and secondary collision", vertical_wall_and_second_hit},
        {"V-groove reflection and reverse path", v_groove_reflection_and_reverse_path},
        {"specular reflection", reflection_geometry}, {"single Cl stoichiometry", chemistry_stoichiometry},
        {"inventory and seed", inventories_and_repeatability}, {"fit cache update", fit_cache_updates},
        {"reaction-cell normal mapping", surface_cell_mapping}, {"periodic cell mapping", periodic_cell_mapping},
        {"curve intersection support", curve_support_is_set_by_the_tracer}
    };
    int failed = 0;
    for (auto& [name, test] : tests) {
        try { test(); std::cout << "PASS " << name << '\n'; }
        catch (const std::exception& e) { ++failed; std::cerr << "FAIL " << name << ": " << e.what() << '\n'; }
    }
    std::cout << tests.size() - failed << "/" << tests.size() << " test groups passed\n";
    return failed ? 1 : 0;
}
