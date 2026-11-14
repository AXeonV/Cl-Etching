#include "planar/backend.hpp"
#include <iostream>
#include <algorithm>
#include <limits>

using namespace etch;
using namespace etch::planar;

// Exhaustive fitted-patch intersections provide a reference independent of the DDA index.
std::optional<Hit> exhaustive_hit(Grid& grid, Surface& surface, Transport& transport, Particle ray) {
    double nearest=std::numeric_limits<double>::infinity();
    std::optional<Hit> best;
    for(int cell:grid.surface_cells) for(int image=-1;image<=1;++image) {
        Curve curve=surface.curve(cell%grid.geometry.nx,cell/grid.geometry.nx);
        if(!curve.defined)continue;
        curve.origin.x += image*grid.geometry.nx;
        double first=1e-8,last=nearest,half=surface.fit_radius()+0.5;
        bool intersects=true;
        for(auto axis:{std::array<double,3>{ray.position.x,ray.direction.x,curve.origin.x},
                       {ray.position.z,ray.direction.z,curve.origin.z}}) {
            auto [p,d,c]=axis;
            if(d==0) { if(p<c-half||p>c+half) intersects=false; }
            else {
                double a=(c-half-p)/d,b=(c+half-p)/d;
                if(a>b)std::swap(a,b);
                first=std::max(first,a);last=std::min(last,b);
            }
        }
        if(!intersects||first>last)continue;
        auto distance=curve.intersection(ray.position,ray.direction,first,last);
        if(!distance||*distance>=nearest)continue;
        Vec point=ray.position+ray.direction*(*distance);
        Vec normal=curve.normal(dot(point-curve.origin,curve.tangent));
        auto owner=transport.reaction_cell(point,normal,curve.origin);
        if(owner){nearest=*distance;best=Hit{point,normal,grid.index(owner->x,owner->z)};}
    }
    return best;
}

int mask_corner_reflections() {
    Grid grid(Geometry{256,512,60,24,72});Surface surface(grid);Transport transport(grid,surface);
    int misses=0,self_hits=0,reference_errors=0;
    for(int i=0;i<80;++i)for(int j=0;j<80;++j) {
        double z=55+(i+0.5)*6/80,dx=0.01+(j+0.5)*0.3/80;
        Vec direction=unit({dx,1});
        Particle ray{{164.5-(z-0.25)*direction.x/direction.z,0.25},direction,100};
        auto first=transport.next_hit(ray);
        if(!first){++misses;continue;}
        ray.position=first->position;ray.direction=reflected(ray.direction,first->normal);
        auto second=transport.next_hit(ray);
        if(!second){if(transport.last_end==RayEnd::MissedSurface)++misses;}
        else if(norm(second->position-first->position)<1e-4)++self_hits;
        if(i%8==0 && j%8==0) {
            auto reference=exhaustive_hit(grid,surface,transport,ray);
            if(bool(reference)!=bool(second) || (reference&&norm(reference->position-second->position)>1e-6))
                ++reference_errors;
        }
    }
    std::cout<<(misses||self_hits||reference_errors?"FAIL ":"PASS ")
             <<"6400 mask-corner reflection rays: misses="<<misses<<", self_hits="<<self_hits
             <<", exhaustive-reference errors="<<reference_errors<<"/100\n";
    return misses||self_hits||reference_errors;
}

int main() {
    // Software diagnostic fixture, independent of the paper's figure inputs.
    Geometry geometry{256, 512, 60, 24, 72};
    Grid flat(geometry); Surface surface(flat); Transport transport(flat, surface);
    int wrong_cells = 0;
    for (int i = 0; i < 100; ++i) {
        double x = 125 + (i + 0.5) / 100;
        auto hit = transport.next_hit({{x, 1}, {0, 1}, 100});
        if (!hit || hit->cell % geometry.nx != static_cast<int>(std::floor(x))) ++wrong_cells;
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
    int corner_errors=mask_corner_reflections();
    return wrong_cells || simulation.transport.interior_exits || corner_errors ? 1 : 0;
}
