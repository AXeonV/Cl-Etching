#include "planar/backend.hpp"
#include "spatial/backend.hpp"
#include <iostream>
#include <algorithm>
#include <vector>
#include <limits>
#include <stdexcept>
using namespace etch;
template<class Backend> void roundtrip(const char* name,typename Backend::Geometry geometry,double dy) {
    typename Backend::Grid grid(geometry);typename Backend::Surface surface(grid);typename Backend::Transport transport(grid,surface);
    int second_count=0,misses=0,missing_back=0,self=0;std::vector<double> error;
    for(int i=0;i<80;++i)for(int j=0;j<80;++j) {
        double target_z=55+(i+0.5)*6/80,dx=0.01+(j+0.5)*0.3/80;
        typename Backend::Source source;
        using Vector=decltype(source.emit(std::declval<Random&>(),geometry).position);
        Vector d,p;
        if constexpr(std::is_same_v<Backend,planar::Backend>) {d=planar::unit({dx,1});p={164.5-(target_z-0.25)*dx,0.25};}
        else {d=spatial::unit({dx,dy,1});p={164.5-(target_z-0.25)*dx,6.25-(target_z-0.25)*dy,0.25};}
        auto first=transport.next_hit({p,d,100});if(!first){++misses;continue;}
        Vector outgoing=reflected(d,first->normal);
        auto second=transport.next_hit({first->position,outgoing,100});
        if(!second) {if(transport.last_end==RayEnd::MissedSurface)++misses;continue;}
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
    if constexpr(std::is_same_v<Backend,planar::Backend>)source={gaussian(0,0.075,-1.5,1.5),constant(100)};
    else source={spatial::gaussian_direction(0.075),constant(100)};
    Random random(42);
    for(int i=0;i<2000;++i)if(!transport.next_hit(source.emit(random,geometry)))throw std::runtime_error("refined reconstruction missed an interface root");
    std::cout<<"PASS refined shared-face interception: 2000 rays\n";
}
int main() {
    refined_interception<planar::Backend>({64,96,20,6,24});
    refined_interception<spatial::Backend>({24,12,64,16,4,12,12});
    cache_check<planar::Backend>({16,32,8,2,8});
    cache_check<spatial::Backend>({16,8,32,8,2,8,8});
    cache_check<planar::Backend>({16,32,8,2,8},2);
    cache_check<spatial::Backend>({16,8,32,8,2,8,8},2);
    roundtrip<planar::Backend>("planar",{256,128,60,24,72},0);
    for(double dy:{0.0,0.1,-0.1})roundtrip<spatial::Backend>("spatial",{256,12,128,60,24,72,12},dy);
}
