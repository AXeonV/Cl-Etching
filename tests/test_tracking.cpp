#include "planar/backend.hpp"
#include "spatial/backend.hpp"
#include <iostream>
#include <algorithm>
#include <vector>
#include <limits>
#include <stdexcept>
#include <type_traits>
using namespace etch;
template<class Backend> void roundtrip(const char* name,typename Backend::Geometry geometry,double dy) {
    typename Backend::Grid grid(geometry);typename Backend::Surface surface(grid);typename Backend::Transport transport(grid,surface);
    int second_count=0,misses=0,missing_back=0,self=0;std::vector<double> error;
    for(int i=0;i<80;++i)for(int j=0;j<80;++j) {
        double target_z=55+(i+0.5)*6/80,dx=0.01+(j+0.5)*0.3/80;
        typename Backend::Source source;
        using Vector=decltype(source.emit(std::declval<Random&>(),geometry).position);
        Vec3 d=unit(Vec3{dx,dy,1});Vector p;
        if constexpr(std::is_same_v<Backend,planar::Backend>) {p={164.5-(target_z-0.25)*dx,0.25};}
        else {p={164.5-(target_z-0.25)*dx,6.25-(target_z-0.25)*dy,0.25};}
        auto first=transport.next_hit({p,d,100});if(!first){++misses;continue;}
        Vec3 outgoing=reflected(d,first->normal);
        auto second=transport.next_hit({first->position,outgoing,100});
        if(!second) {if(transport.diagnostics.last_end==RayEnd::MissedSurface)++misses;continue;}
        ++second_count;
        if(norm(second->position-first->position)<1e-4)++self;
        auto back=transport.next_hit({second->position,outgoing*-1,100});
        if(!back){++missing_back;error.push_back(std::numeric_limits<double>::infinity());continue;}
        Vector delta=back->position-first->position;delta.x-=geometry.nx*std::round(delta.x/geometry.nx);
        if constexpr(std::is_same_v<Backend,spatial::Backend>)delta.y-=geometry.ny*std::round(delta.y/geometry.ny);
        error.push_back(norm(delta));
    }
    std::sort(error.begin(),error.end());
    if(misses||missing_back||self||error.back()>1e-7)throw std::runtime_error("reverse trajectory consistency failed");
    std::cout<<name<<" dy="<<dy<<" rays=6400 second="<<second_count<<" misses="<<misses<<" immediate_self="<<self<<" missing_back="<<missing_back;
    for(double threshold:{1e-5,0.5,1.0,2.0})std::cout<<" error_above_"<<threshold<<'='<<std::count_if(error.begin(),error.end(),[&](double e){return e>threshold;});
    std::cout<<" p50="<<error[error.size()/2]<<" p90="<<error[error.size()*9/10]<<" p99="<<error[error.size()*99/100]<<'\n';
}
template<class Backend> void cache_check(typename Backend::Geometry geometry,int n=1) {
    typename Backend::Grid grid(geometry);typename Backend::Surface surface(grid,3,n);
    for(int z=0;z<geometry.nz*n;++z)for(int x=0;x<geometry.nx*n;++x) {
        if constexpr(std::is_same_v<Backend,planar::Backend>)surface.field(x,z);
        else for(int y=0;y<geometry.ny*n;++y)surface.field(x,y,z);
    }
    if constexpr(std::is_same_v<Backend,planar::Backend>){grid.remove_empty_cell(8,8);surface.invalidate(8,8);}
    else {grid.remove_empty_cell(grid.index(8,4,8));surface.invalidate(8,4,8);}
    typename Backend::Surface fresh(grid,3,n);double largest=0;
    for(int z=0;z<geometry.nz*n;++z)for(int x=0;x<geometry.nx*n;++x) {
        if constexpr(std::is_same_v<Backend,planar::Backend>) {
            largest=std::max(largest,std::abs(surface.field(x,z).value({0.37,0.61})-fresh.field(x,z).value({0.37,0.61})));
        }else for(int y=0;y<geometry.ny*n;++y)largest=std::max(largest,std::abs(surface.field(x,y,z).value({0.37,0.23,0.61})-fresh.field(x,y,z).value({0.37,0.23,0.61})));
    }
    if(largest>1e-12)throw std::runtime_error("stale field after cell update");
    std::cout<<"PASS field cache vs fresh reconstruction: "<<largest<<'\n';
}
template<class Backend> void refined_interception(typename Backend::Geometry geometry) {
    typename Backend::Grid grid(geometry);typename Backend::Surface surface(grid,3,4);typename Backend::Transport transport(grid,surface);
    typename Backend::Source source;
    source={gaussian_direction(0.075),constant(100)};
    Random random(42);
    for(int i=0;i<2000;++i)if(!transport.next_hit(source.emit(random,geometry)))throw std::runtime_error("refined reconstruction missed an interface root");
    std::cout<<"PASS refined shared-face interception: 2000 rays\n";
}
void extruded_2d3v_trajectories() {
    static_assert(std::is_same_v<decltype(planar::Particle{}.position),planar::Vec>);
    static_assert(std::is_same_v<decltype(planar::Particle{}.direction),Vec3>);
    static_assert(std::is_same_v<decltype(planar::Hit{}.normal),Vec3>);
    planar::Grid plane({40,80,20,6,16});
    spatial::Grid volume({40,12,80,20,6,16,12});
    for(int z=20;z<50;++z)for(int x=12;x<28;++x) {
        int cell=plane.index(x,z);plane.material[cell]=Gas;plane.sites[cell]={};
        for(int y=0;y<12;++y) {
            cell=volume.index(x,y,z);volume.material[cell]=Gas;volume.sites[cell]={};
        }
    }
    plane.refresh(20,40,80);volume.refresh(20,6,40,80);
    planar::Surface curve(plane);spatial::Surface face(volume);
    planar::Transport planar_transport(plane,curve);spatial::Transport spatial_transport(volume,face);
    int compared=0;
    for(double dy:{-4.0,-0.3,0.0,0.5,3.0}) {
        Vec3 direction=unit(Vec3{1,dy,0.2});
        planar::Particle a{{20,30},direction,100};
        spatial::Particle b{{20,6.25,30},direction,100};
        for(int bounce=0;bounce<4;++bounce) {
            auto first=planar_transport.next_hit(a);auto second=spatial_transport.next_hit(b);
            if(!first||!second)throw std::runtime_error("extruded trajectories must hit both walls");
            if(norm(first->position-planar::projection(second->position))>1e-8 ||
               norm(first->normal-second->normal)>1e-8 || first->normal.y!=0)
                throw std::runtime_error("2D3V and extruded 3D collision geometry disagree");
            double mu_a=-dot(a.direction,first->normal),mu_b=-dot(b.direction,second->normal);
            if(std::abs(mu_a-mu_b)>1e-8)throw std::runtime_error("three-component incidence angle mismatch");
            a.position=first->position;b.position=second->position;
            plane.wrap_position(a.position);volume.wrap_position(b.position);
            a.direction=reflected(a.direction,first->normal);b.direction=reflected(b.direction,second->normal);
            if(norm(a.direction-b.direction)>1e-8||std::abs(norm(a.direction)-1)>1e-12||a.direction.y!=direction.y)
                throw std::runtime_error("2D3V reflection does not preserve the three-component direction");
            ++compared;
        }
    }
    std::cout<<"PASS 2D3V vs extruded 3D: "<<compared<<" wall collisions and reflections\n";
}
int main() {
    extruded_2d3v_trajectories();
    refined_interception<planar::Backend>({64,96,20,6,24});
    refined_interception<spatial::Backend>({24,12,64,16,4,12,12});
    cache_check<planar::Backend>({16,32,8,2,8});
    cache_check<spatial::Backend>({16,8,32,8,2,8,8});
    cache_check<planar::Backend>({16,32,8,2,8},2);
    cache_check<spatial::Backend>({16,8,32,8,2,8,8},2);
    for(double dy:{0.0,0.3,-2.0})roundtrip<planar::Backend>("planar",{256,128,60,24,72},dy);
    for(double dy:{0.0,0.1,-0.1})roundtrip<spatial::Backend>("spatial",{256,12,128,60,24,72,12},dy);
}
