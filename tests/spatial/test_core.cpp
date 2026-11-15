#include "spatial/backend.hpp"
#include "spatial/traversal.hpp"
#include "common/linear.hpp"
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <set>
using namespace etch;
using namespace etch::spatial;
void check(bool value,const std::string& message) { if(!value)throw std::runtime_error(message); }
void close(double actual,double expected,double tolerance=1e-6) {
    check(std::abs(actual-expected)<tolerance,"expected "+std::to_string(expected)+", got "+std::to_string(actual));
}
Geometry geometry() { return {32,32,80,20,6,16,32}; }
void vectors() {
    Vec x{1,0,0},y{0,1,0};close(norm(cross(x,y)-Vec{0,0,1}),0);
    Vec d=unit({1,2,3}),n=unit({2,-1,1}),r=reflected(d,n);
    close(norm(r),1);close(dot(r,n),-dot(d,n));close(norm(reflected(r,n)-d),0);
}
void voxel_storage_and_hole() {
    Geometry g=geometry();g.opening_y=12;Grid grid(g);
    check(grid.index(-1,-1,20)==grid.index(31,31,20),"both lateral axes are periodic");
    check(grid.coordinates(grid.index(5,7,30)).y==7,"z-y-x storage coordinates");
    check(grid.at(16,16,18)==Gas && grid.at(16,2,18)==Mask,"finite rectangular mask aperture");
    check(grid.at(16,16,20)==Silicon && grid.boundary[grid.index(16,16,20)],"six-face surface cells");
    check(grid.silicon_inventory()==32LL*32*60*30,"initial three-dimensional inventory");
}
void supercover_vertices_and_edges() {
    VoxelWalk walk({0.5,0.5,0.5},unit({1,1,1}));std::set<std::array<int,3>> visited;
    for(int i=0;i<7;++i){walk.advance();visited.insert(walk.cell);}
    check(visited.size()==7 && visited.count({1,1,1}),"visit all seven vertex-touching cells");
    for(int mask=1;mask<8;++mask)check(visited.count({mask&1,(mask>>1)&1,(mask>>2)&1}),"vertex supercover subset");
    VoxelWalk edge({0.5,0.5,0.5},unit({1,1,0}));visited.clear();
    for(int i=0;i<3;++i){edge.advance();visited.insert(edge.cell);}
    check(visited.size()==3 && visited.count({1,1,0}),"edge crossing visits three cells");
}
Patch flat() { Patch p;p.origin={0,0,10};p.u_axis={1,0,0};p.v_axis={0,1,0};p.outward={0,0,-1};return p; }
void quadratic_patch() {
    Patch p=flat();p.a=0.1;p.b=0.2;p.c=0.05;p.d=0.3;p.e=0.1;p.f=0.4;
    close(p.height(1,2),1.9);close(norm(p.position(1,2)-Vec{1,2,8.1}),0);
    Vec expected=unit(Vec{-0.6,-0.95,-1});close(norm(p.normal(1,2)-expected),0);
}
void normal_projection() {
    Patch p=flat();auto projected=p.project({2,3,7});close(projected.u,2);close(projected.v,3);close(projected.distance,3);
    p=flat();p.a=1;p.b=1;p.origin={0,0,0};p.outward={0,0,1};
    auto q=p.project({0,0,2});close(q.u*q.u+q.v*q.v,1.5,1e-5);close(q.distance,std::sqrt(1.75),1e-5);
    Vec delta=Vec{0,0,2}-p.position(q.u,q.v);
    close(dot(delta,Vec{1,0,2*q.u}),0,1e-5);close(dot(delta,Vec{0,1,2*q.v}),0,1e-5);
}
void six_coefficient_fit() {
    double matrix[6][6]={},rhs[6]={},result[6]={},expected[]={0.1,-0.2,0.05,0.4,-0.3,0.7};
    for(int u=-3;u<=3;++u)for(int v=-3;v<=3;++v){
        double row[]={double(u*u),double(v*v),double(u*v),double(u),double(v),1},w=0;
        for(int i=0;i<6;++i)w+=row[i]*expected[i];
        for(int i=0;i<6;++i){rhs[i]+=row[i]*w;for(int j=0;j<6;++j)matrix[i][j]+=row[i]*row[j];}
    }
    check(least_squares(matrix,rhs,result)==6,"full rank six-coefficient fit");
    for(int i=0;i<6;++i)close(result[i],expected[i],1e-9);
}
void fitted_plane_and_periodic_mapping() {
    Geometry g=geometry();g.opening=g.nx;Grid grid(g);Surface surface(grid);Transport transport(grid,surface);
    const auto& patch=surface.patch(16,16,20);check(patch.defined,"flat patch is defined");close(norm(patch.normal(0,0)-Vec{0,0,-1}),0);
    auto hit=transport.next_hit({{16.3,-0.2,1},{0,0,1},100});
    check(hit.has_value(),"ray crossing y seam hits plane");close(hit->position.z,20.5);
    Cell owner=grid.coordinates(hit->cell);check(owner.x==16&&owner.y==31&&owner.z==20,"periodic y reaction cell");
    auto mapped=transport.reaction_cell({16.3,16.7,19.5},{0,0,-1},{16.5,16.5,20.5});
    check(mapped&&mapped->z==20,"normal maps gas point to surface voxel");
    check(!transport.next_hit({{16,16,1},{0,1,0},100})&&transport.last_end==RayEnd::PeriodicOrbit,"lateral periodic orbit");
    check(!transport.next_hit({{16,16,19},unit({1,2,0}),100})&&transport.last_end==RayEnd::PeriodicOrbit,"arbitrary horizontal ray above a plane terminates");
}
void closed_rational_orbit() {
    Geometry geometry{32,32,64,30,0,32,32};Grid grid(geometry);
    for(int z=12;z<19;++z)for(int y=15;y<17;++y)for(int x=15;x<17;++x){int cell=grid.index(x,y,z);grid.material[cell]=Silicon;grid.sites[cell][0]=30;}
    grid.refresh(16,16,32,64);Surface surface(grid);Transport transport(grid,surface);
    auto hit=transport.next_hit({{0.25,0.25,15.25},unit({1,2,0}),100});
    check(!hit&&transport.last_end==RayEnd::PeriodicOrbit,"clear rational periodic ray terminates without an arbitrary bounce cap");
}
void reflection_in_y_plane() {
    Geometry g=geometry();g.opening=g.nx;Grid grid(g);
    for(int z=20;z<50;++z)for(int y=6;y<22;++y)for(int x=0;x<32;++x){
        int cell=grid.index(x,y,z);grid.material[cell]=Gas;grid.sites[cell]={};
    }
    grid.refresh(16,16,40,80);Surface surface(grid);Transport transport(grid,surface);
    Particle ray{{16,20,30},unit({0,1,2}),100};auto wall=transport.next_hit(ray);
    check(wall.has_value(),"y-facing wall hit");close(wall->position.y,22.5);close(wall->position.z,35);
    close(norm(wall->normal-Vec{0,-1,0}),0);
    ray.position=wall->position;ray.direction=reflected(ray.direction,wall->normal);
    auto floor=transport.next_hit(ray);check(floor.has_value(),"secondary floor hit after y reflection");
    close(floor->position.y,14.75);close(floor->position.z,50.5);
    auto reverse=transport.next_hit({floor->position,ray.direction*-1,100});
    check(reverse.has_value(),"reverse three-dimensional trajectory");close(norm(reverse->position-wall->position),0);
}
void corner_normal_mapping() {
    Geometry g{32,32,64,0,0,32,32};Grid grid(g);
    for(int z=0;z<64;++z)for(int y=0;y<32;++y)for(int x=0;x<32;++x)if(x+y+z<75){
        int cell=grid.index(x,y,z);grid.material[cell]=Gas;grid.sites[cell]={};
    }
    grid.refresh(16,16,32,64);Surface surface(grid);Transport transport(grid,surface);
    auto mapped=transport.reaction_cell({25,25,26.5},unit({-1,-1,-1}),{25,25,26.5});
    check(mapped&&grid.boundary[grid.index(mapped->x,mapped->y,mapped->z)],"normal crossing an edge finds touching surface cell");
}
void neighbor_transfer_and_cache() {
    for(int y:{16,0}) {
        Geometry g=geometry();g.opening=g.nx;Grid grid(g);int original=grid.index(16,y,20),tip=grid.index(16,y-1,19);
        grid.sites[original]={1,1,1,1};grid.material[tip]=Silicon;Sites inventory{10,4,5,6};grid.sites[tip]=inventory;
        grid.refresh(16,y,20,3);Surface surface(grid);CellUpdater updater(grid,surface);
        surface.candidates(16,y,19);surface.patch(16,y,20);auto fits=surface.fit_count;
        long long si=grid.silicon_inventory(),cl=grid.chlorine_inventory();
        auto result=updater.remove(original);
        check(result.relocated&&grid.material[tip]==Gas&&grid.sites[original]==inventory,"26-neighbor transfer in y, including periodic seam");
        check(grid.silicon_inventory()+4==si&&grid.chlorine_inventory()+6==cl,"preserve Si and Cl during relocation");
        surface.patch(16,y,20);check(surface.fit_count>fits,"invalidate affected bivariate fit");
        const auto& cells=surface.candidates(16,y,19);
        check(std::find(cells.begin(),cells.end(),tip)==cells.end(),"remove deleted voxel from candidate index");
    }
}
void source_solid_angle_distributions() {
    Random rng(42);auto cosine=cosine_direction(),gaussian=gaussian_direction(0.07);
    double mu=0,mu2=0,x2=0,y2=0,theta2=0;
    for(int i=0;i<80000;++i){Vec d=cosine(rng);close(norm(d),1);mu+=d.z;mu2+=d.z*d.z;x2+=d.x*d.x;y2+=d.y*d.y;
        Vec g=gaussian(rng);check(g.z>0,"downward Gaussian hemisphere");theta2+=std::pow(std::acos(g.z),2);}
    close(mu/80000,2.0/3,0.004);close(mu2/80000,0.5,0.004);close(x2/80000,0.25,0.004);close(y2/80000,0.25,0.004);
    close(theta2/80000,2*0.07*0.07,0.0003);
}
void shared_engine_inventory() {
    Geometry g{24,12,64,16,4,12,12};Settings s{g};s.max_bounces=500;
    for(int model=0;model<2;++model) {
        auto make=[&]()->std::unique_ptr<Reaction>{if(model==0)return std::make_unique<YieldReaction>();return std::make_unique<ChlorinationReaction>();};
        Simulation a(s,Source{gaussian_direction(0.04),constant(100)},make());
        Simulation b(s,Source{gaussian_direction(0.04),constant(100)},make());
        a.advance(3000);b.advance(3000);
        check(a.grid.sites==b.grid.sites&&a.grid.material==b.grid.material,"repeatable spatial trajectories");
        check(a.grid.silicon_inventory()+a.reaction->removed+a.stats.pruned_si==a.initial_si,"spatial Si inventory");
        if(model==1)check(a.reaction->captures==a.grid.chlorine_inventory()+4*a.reaction->products_sicl4+a.stats.pruned_cl,"spatial Cl inventory");
        check(a.transport.interior_exits==0,"spatial beam does not pass through substrate");
    }
}
void centroid_degeneracy() {
    Geometry g=geometry();Grid grid(g);int cell=grid.index(16,16,10);grid.material[cell]=Silicon;grid.sites[cell][0]=12;
    grid.refresh(16,16,10,2);Surface surface(grid);
    check(!surface.patch(16,16,10).defined&&surface.undefined_normal_fits==1,"isolated voxel has no invented normal");
    check(grid.remaining(cell)==12,"retain isolated voxel inventory");
}
int main() {
    std::cout<<std::unitbuf;
    std::vector<std::pair<std::string,std::function<void()>>> tests={
        {"spatial vectors",vectors},{"voxel storage and aperture",voxel_storage_and_hole},
        {"vertex and edge supercover",supercover_vertices_and_edges},{"local quadratic surface geometry",quadratic_patch},
        {"normal projection",normal_projection},{"six coefficient least squares",six_coefficient_fit},
        {"fitted plane and y periodicity",fitted_plane_and_periodic_mapping},{"rational periodic orbit",closed_rational_orbit},{"y reflection and reverse path",reflection_in_y_plane},
        {"corner normal mapping",corner_normal_mapping},{"26-neighbor inventory and cache",neighbor_transfer_and_cache},
        {"solid-angle source distributions",source_solid_angle_distributions},{"shared engine inventory",shared_engine_inventory},
        {"centroid degeneracy",centroid_degeneracy}
    };
    int failed=0;for(auto& [name,test]:tests){try{test();std::cout<<"PASS "<<name<<'\n';}
        catch(const std::exception& error){++failed;std::cerr<<"FAIL "<<name<<": "<<error.what()<<'\n';}}
    std::cout<<tests.size()-failed<<'/'<<tests.size()<<" spatial test groups passed\n";return failed?1:0;
}
