#include <cfenv>
#include "spatial/backend.hpp"
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
CellRemoval update_cell(CellUpdater& updater,const Grid& grid,int cell,Vec3 normal) {
    auto p=grid.coordinates(cell);return updater.remove(cell,normal,{p.x+0.5,p.y+0.5,p.z+0.5});
}
void voxel_storage_and_hole() {
    Geometry g=geometry();g.opening_y=12;Grid grid(g);
    check(grid.index(-1,-1,20)==grid.index(31,31,20),"both lateral axes are periodic");
    check(grid.coordinates(grid.index(5,7,30)).y==7,"z-y-x storage coordinates");
    check(grid.at(16,16,18)==Gas && grid.at(16,2,18)==Mask,"finite rectangular mask aperture");
    check(grid.at(16,16,20)==Silicon && grid.boundary[grid.index(16,16,20)],"six-face surface cells");
    check(grid.silicon_inventory()==32LL*32*60*30,"initial three-dimensional inventory");
}
Patch flat() { Patch p;p.origin={0,0,10};p.u_axis={1,0,0};p.v_axis={0,1,0};p.outward={0,0,-1};return p; }
void quadratic_patch() {
    Patch p=flat();p.a=0.1;p.b=0.2;p.c=0.05;p.d=0.3;p.e=0.1;p.f=0.4;
    close(p.height(1,2),1.9);close(norm(p.position(1,2)-Vec3{1,2,8.1}),0);
    Vec3 expected=unit(Vec3{-0.6,-0.95,-1});close(norm(p.normal(1,2)-expected),0);
}
void normal_projection() {
    Patch p=flat();auto projected=p.project({2,3,7});close(projected.u,2);close(projected.v,3);close(projected.distance,3);
    p=flat();p.a=1;p.b=1;p.origin={0,0,0};p.outward={0,0,1};
    auto q=p.project({0,0,2});close(q.u*q.u+q.v*q.v,1.5,1e-5);close(q.distance,std::sqrt(1.75),1e-5);
    Vec3 delta=Vec3{0,0,2}-p.position(q.u,q.v);
    close(dot(delta,Vec3{1,0,2*q.u}),0,1e-5);close(dot(delta,Vec3{0,1,2*q.v}),0,1e-5);
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
    const auto& patch=surface.patch(16,16,20);check(patch.defined,"flat patch is defined");close(norm(patch.normal(0,0)-Vec3{0,0,-1}),0);
    auto hit=transport.next_hit({{16.3,-0.2,1},{0,0,1},100});
    check(hit.has_value(),"ray crossing y seam hits plane");close(hit->position.z,20.5);
    Cell owner=grid.coordinates(hit->cell);check(owner.x==16&&owner.y==31&&owner.z==20,"periodic y reaction cell");
    auto mapped=transport.reaction_cell({16.3,16.7,19.5},{0,0,-1},{16.5,16.5,20.5});
    check(mapped&&mapped->z==20,"normal maps gas point to surface voxel");
    check(!transport.next_hit({{16,16,1},{0,1,0},100})&&transport.diagnostics.last_end==RayEnd::PeriodicOrbit,"lateral periodic orbit");
    check(!transport.next_hit({{16,16,19},unit({1,2,0}),100})&&transport.diagnostics.last_end==RayEnd::PeriodicOrbit,"arbitrary horizontal ray above a plane terminates");
}
void closed_rational_orbit() {
    Geometry geometry{32,32,64,30,0,32,32};Grid grid(geometry);
    for(int z=12;z<19;++z)for(int y=15;y<17;++y)for(int x=15;x<17;++x){int cell=grid.index(x,y,z);grid.material[cell]=Silicon;grid.sites[cell][0]=30;}
    grid.refresh(16,16,32,64);Surface surface(grid);Transport transport(grid,surface);
    auto hit=transport.next_hit({{0.25,0.25,15.25},unit({1,2,0}),100});
    check(!hit&&transport.diagnostics.last_end==RayEnd::PeriodicOrbit,"clear rational periodic ray terminates without an arbitrary bounce cap");
}
void reflection_in_y_plane() {
    Geometry g=geometry();g.opening=g.nx;Grid grid(g);
    for(int z=20;z<50;++z)for(int y=6;y<22;++y)for(int x=0;x<32;++x){
        int cell=grid.index(x,y,z);grid.material[cell]=Gas;grid.sites[cell]={};
    }
    grid.refresh(16,16,40,80);Surface surface(grid);Transport transport(grid,surface);
    Particle ray{{16,20,30},unit({0,1,2}),100};auto wall=transport.next_hit(ray);
    check(wall.has_value(),"y-facing wall hit");close(wall->position.y,22.5);close(wall->position.z,35);
    close(norm(wall->normal-Vec3{0,-1,0}),0);
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
        int intermediate=grid.index(16,y-1,20);Sites near_inventory{6,7,8,9};grid.sites[intermediate]=near_inventory;
        grid.refresh(16,y,20,3);Surface surface(grid);CellUpdater updater(grid,surface);
        surface.candidates(16,y,19);surface.patch(16,y,20);auto fits=surface.fit_count;
        long long si=grid.silicon_inventory(),cl=grid.chlorine_inventory();
        auto result=update_cell(updater,grid,original,{0,0,-1});
        check(result.relocated&&grid.material[tip]==Gas&&grid.sites[original]==Sites{11,5,6,7}&&grid.sites[intermediate]==near_inventory,"connected donor transfer in y, including periodic seam");
        check(grid.silicon_inventory()==si&&grid.chlorine_inventory()==cl,"preserve Si and Cl during relocation");
        surface.patch(16,y,20);check(surface.fit_count>fits,"invalidate affected bivariate fit");
        const auto& cells=surface.candidates(16,y,19);
        check(std::find(cells.begin(),cells.end(),tip)==cells.end(),"remove deleted voxel from candidate index");
    }
}
void source_solid_angle_distributions() {
    Random rng(42);auto cosine=cosine_direction(),gaussian=gaussian_direction(0.07);
    double mu=0,mu2=0,x2=0,y2=0,theta2=0;
    for(int i=0;i<80000;++i){Vec3 d=cosine(rng);close(norm(d),1);mu+=d.z;mu2+=d.z*d.z;x2+=d.x*d.x;y2+=d.y*d.y;
        Vec3 g=gaussian(rng);check(g.z>0,"downward Gaussian hemisphere");theta2+=std::pow(std::acos(g.z),2);}
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
        check(a.grid.silicon_inventory()+a.reaction->removed+a.stats.detached_si==a.initial_si,"spatial Si inventory");
        if(model==1)check(a.reaction->captures==a.grid.chlorine_inventory()+4*a.reaction->products_sicl4+a.stats.detached_cl,"spatial Cl inventory");
        check(a.transport.diagnostics.interior_exits==0,"spatial beam does not pass through substrate");
    }
}
void centroid_degeneracy() {
    Geometry g=geometry();Grid grid(g);int cell=grid.index(16,16,10);grid.material[cell]=Silicon;grid.sites[cell][0]=12;
    grid.refresh(16,16,10,2);Surface surface(grid);
    check(!surface.patch(16,16,10).defined&&surface.undefined_normal_fits==1,"isolated voxel has no invented normal");
    check(grid.remaining(cell)==12,"retain isolated voxel inventory");
}
void connected_donor_selection() {
    Random selection_rng(17);
    for(int axis:{0,1})for(int coordinate:{8,0}) {
        Grid grid({16,16,24,20,0,16,16});
        std::fill(grid.material.begin(),grid.material.end(),Gas);
        auto cell=[&](int offset,int z){return axis==0?grid.index(coordinate+offset,8,z):grid.index(8,coordinate+offset,z);};
        int bridge=cell(0,8),leaf=cell(-1,8);
        for(int offset:{-1,0,1})grid.material[cell(offset,8)]=Silicon;
        RemovalConnectivity<Grid> connectivity(grid);
        check(!connectivity.preserves_connection(bridge),"a three-dimensional neck is not a removable leaf");
        std::vector<RemovalCandidate> candidates{{bridge,-1},{leaf,0}};
        check(connectivity.select(candidates,selection_rng).cell==leaf&&connectivity.select(candidates,selection_rng).redirected,"prefer an available deletion that preserves connectivity");
        candidates={{bridge,0}};
        check(connectivity.select(candidates,selection_rng).cell==bridge&&!connectivity.select(candidates,selection_rng).connected,"unavoidable fracture is not pinned");
        for(int z=4;z<8;++z)for(int offset:{-1,1})grid.material[cell(offset,z)]=Silicon;
        for(int offset=-1;offset<=1;++offset)grid.material[cell(offset,4)]=Silicon;
        check(connectivity.preserves_connection(bridge),"global alternate path is respected across either periodic seam");
    }
    Grid grid({16,16,24,20,0,16,16});std::fill(grid.material.begin(),grid.material.end(),Gas);
    int original=grid.index(8,8,8),attached=grid.index(7,8,8),diagonal=grid.index(9,9,8),mask=grid.index(9,8,8);
    grid.material[original]=grid.material[attached]=grid.material[diagonal]=Silicon;grid.material[mask]=Mask;
    RemovalConnectivity<Grid> connectivity(grid);
    auto branch=connectivity.local_silicon_branch(original,{original,attached,diagonal,mask});
    check(branch.cells.size()==2&&branch.contains(attached),"six-face inventory branch excludes diagonal material across mask");
}

void collision_normal_at_symmetric_branch() {
    for(int sign:{-1,1}) {
        Grid grid({20,20,32,24,0,20,20});
        std::fill(grid.material.begin(),grid.material.end(),Gas);
        std::fill(grid.sites.begin(),grid.sites.end(),Sites{});
        auto cell=[&](int x){return grid.index(x,8,8);};
        for(int x=7;x<=9;++x){grid.material[cell(x)]=Silicon;grid.sites[cell(x)]={26,0,0,0};}
        grid.sites[cell(8)]={1,1,1,1};grid.refresh(8,8,8,4);
        Surface surface(grid);CellUpdater updater(grid,surface);
        auto result=update_cell(updater,grid,cell(8),{double(sign),0,0});
        check(result.cell==cell(8+sign),"collision normal orients removal when the local centroid is symmetric");
        check(updater.diagnostics.undefined_fit_removals==1&&updater.diagnostics.unavoidable_splits==0,"singular local fit retains surface connection");
    }
}

void fragment_release_and_anchors() {
    Grid grid({20,20,32,24,0,20,20});
    std::fill(grid.material.begin(),grid.material.end(),Gas);
    std::fill(grid.sites.begin(),grid.sites.end(),Sites{});
    auto cell=[&](int z){return grid.index(8,8,z);};
    for(int z=5;z<=11;++z){grid.material[cell(z)]=Silicon;grid.sites[cell(z)]={6,7,8,9};}
    grid.material[cell(12)]=Mask;grid.sites[cell(8)]={};grid.refresh(8,8,8,8);
    Surface surface(grid);CellUpdater updater(grid,surface);
    for(int zz=4;zz<=9;++zz)for(int yy=7;yy<=9;++yy)for(int xx=7;xx<=9;++xx)surface.node(xx,yy,zz);
    auto si=grid.silicon_inventory(),cl=grid.chlorine_inventory();
    auto removal=update_cell(updater,grid,cell(8),{0,0,-1});
    check(removal.cell==cell(7)&&removal.detached_cells==2,"a real neck break releases the unanchored component");
    check(removal.detached_si==60&&removal.detached_cl==100,"released fragment states are accounted separately from reaction products");
    check(grid.silicon_inventory()+removal.detached_si==si,"Si closes after fracture");
    check(grid.chlorine_inventory()+removal.detached_cl==cl,"bound Cl closes after fracture");
    check(grid.material[cell(5)]==Gas&&grid.material[cell(6)]==Gas&&grid.material[cell(12)]==Mask,"release only the unsupported branch");
    Surface fresh(grid);
    for(int zz=4;zz<=9;++zz)for(int yy=7;yy<=9;++yy)for(int xx=7;xx<=9;++xx) {
        auto cached=surface.node(xx,yy,zz),expected=fresh.node(xx,yy,zz);
        for(size_t derivative=0;derivative<cached.size();++derivative)close(cached[derivative],expected[derivative],1e-12);
    }
    RemovalConnectivity<Grid> connectivity(grid);
    check(connectivity.detached_components(cell(7)).empty(),"fragment release reaches a fixed point");
    grid.material[cell(5)]=Mask;
    for(int z:{6,7})grid.material[cell(z)]=Silicon;
    grid.material[cell(8)]=Gas;
    check(connectivity.detached_components(cell(8)).empty(),"separate components attached to different fixed anchors both survive");
}

void conservative_update_defers_geometry() {
    Grid grid({20,20,32,12,0,20,20});int original=grid.index(10,10,12);grid.sites[original]={1,1,1,1};
    Surface surface(grid);CellUpdater updater(grid,surface,0.3,17);
    surface.patch(10,10,12);auto fits=surface.fit_count;auto before=grid.sites;auto material=grid.material;
    auto pending=update_cell(updater,grid,original,{0,0,-1});
    check(pending.reaction_cells==0&&pending.detached_cells==0,"capacity shortage causes no geometric or inventory loss");
    check(grid.sites==before&&grid.material==material&&surface.fit_count==fits,"deferred updates preserve inventory and geometry caches");
    check(updater.diagnostics.deferred_updates==1,"defer is recorded separately from a committed deletion");
    grid.sites[original]={};
    auto committed=update_cell(updater,grid,original,{0,0,-1});
    check(committed.reaction_cells==1&&grid.material[original]==Gas,"empty cell still advances after a deferred update");
}

void unbiased_candidate_ties_and_capacity_support() {
    Grid grid({20,20,32,12,0,20,20});
    auto cell=[&](int x){return grid.index(x,10,12);};
    RemovalConnectivity<Grid> connectivity(grid);Random rng(23);
    int first=0;
    for(int i=0;i<1000;++i) {
        std::vector<RemovalCandidate> candidates{{cell(9),0,1},{cell(9),0,1},{cell(11),0,1}};
        first+=connectivity.select(candidates,rng).cell==cell(9);
    }
    check(first>400&&first<600,"periodic aliases do not bias otherwise equivalent deletion choices");
    std::vector<RemovalCandidate> candidates{{cell(9),0,4},{cell(11),0,1}};auto before_rng=rng;
    check(connectivity.select(candidates,rng).cell==cell(11)&&rng==before_rng,"minimize transport among equal geometric scores without unnecessary RNG draws");
    grid.sites[cell(10)]={1,1,1,1};grid.sites[cell(12)]={25,0,0,0};
    Surface surface(grid);
    CellUpdater local(grid,surface,0.3,17);
    auto pending=update_cell(local,grid,cell(10),{0,0,-1});
    check(pending.reaction_cells==0,"strictly local remapping defers when its own volume is full");

}

void connectivity_matches_global_components() {
    Grid grid({5,4,7,4,0,5,4});Random rng(77);RemovalConnectivity<Grid> connectivity(grid);
    auto components=[&](int removed) {
        std::vector<bool> seen(grid.material.size());int count=0;
        for(size_t seed=0;seed<grid.material.size();++seed)if(int(seed)!=removed&&grid.material[seed]!=Gas&&!seen[seed]) {
            ++count;std::vector<int> queue{int(seed)};seen[seed]=true;
            for(size_t i=0;i<queue.size();++i)for(int neighbor:grid.face_neighbors(queue[i])) {
                if(neighbor>=0&&neighbor!=removed&&grid.material[neighbor]!=Gas&&!seen[neighbor]) {seen[neighbor]=true;queue.push_back(neighbor);}
            }
        }
        return count;
    };
    for(int sample=0;sample<500;++sample) {
        for(auto& material:grid.material){double value=uniform(rng);material=value<0.4?Silicon:value<0.45?Mask:Gas;}
        int before=components(-1),checked=0;
        for(size_t n=0;n<grid.material.size()&&checked<5;++n) {
            int cell=int((n+sample)%grid.material.size());if(grid.material[cell]!=Silicon)continue;
            check(connectivity.preserves_connection(cell)==(components(cell)<=before),"multi-source decision equals full component enumeration");
            ++checked;
        }
    }
}

void collision_projection_and_periodic_image() {
    std::array<int,2> chosen{};
    for(int side=0;side<2;++side) {
        Grid grid({20,20,80,10,0,20,20});
        for(int z=0;z<80;++z)for(int y=0;y<20;++y)for(int x=0;x<20;++x) {
            int cell=grid.index(x,y,z),top=std::max(0,32-(x-10)*(x-10));
            grid.material[cell]=z>=top?Silicon:Gas;grid.sites[cell]=z>=top?Sites{30,0,0,0}:Sites{};
        }
        for(int y=0;y<20;++y)for(int x:{9,11}){int cell=grid.index(x,y,30);grid.material[cell]=Silicon;grid.sites[cell]={30,0,0,0};}
        int original=grid.index(10,10,32);grid.sites[original]={};grid.refresh(10,10,40,80);
        Surface surface(grid);Patch patch=surface.patch(10,10,32);Vec3 offset{side==0?-0.3:0.3,0,0};
        Vec3 point=patch.position(dot(offset,patch.u_axis),dot(offset,patch.v_axis));
        CellUpdater updater(grid,surface);chosen[side]=updater.remove(original,{0,0,-1},point).cell;
    }
    check(chosen[0]!=chosen[1]&&chosen[0]%20+chosen[1]%20==20,"three-dimensional contact projections select reflected counterparts");
    Geometry g{20,20,60,20,0,20,20};Grid a(g),b(g);
    for(auto* grid:{&a,&b}) {
        int tip=grid->index(-1,-1,19);grid->material[tip]=Silicon;grid->sites[tip]={30,0,0,0};
        grid->sites[grid->index(0,0,20)]={};grid->refresh(0,0,20,4);
    }
    Surface sa(a),sb(b);CellUpdater ua(a,sa),ub(b,sb);
    ua.remove(a.index(0,0,20),{0,0,-1},{0.25,0.25,20.5});
    ub.remove(b.index(0,0,20),{0,0,-1},{20.25,20.25,20.5});
    check(a.material==b.material&&a.sites==b.sites,"both periodic contact images give identical updates");
}

void undefined_normal_removal() {
    for (bool symmetric_neighbors : {false, true}) {
        Grid grid({20,16,40,15,4,10,16});
        int x=10,y=8,z=6,original=grid.index(x,y,z);
        grid.material[original]=Silicon;grid.sites[original]={1,1,1,1};
        if(symmetric_neighbors)for(int dx:{-1,1}) {
            int cell=grid.index(x+dx,y,z);grid.material[cell]=Silicon;grid.sites[cell]={8,7,6,5};
        }
        grid.refresh(x,y,z,3);Surface surface(grid);CellUpdater updater(grid,surface);
        check(!surface.patch(x,y,z).defined,"symmetric neighborhood has no centroid normal");
        auto si=grid.silicon_inventory(),cl=grid.chlorine_inventory();
        std::feclearexcept(FE_ALL_EXCEPT);
        auto result=update_cell(updater,grid,original,{0,0,-1});
        check(!std::fetestexcept(FE_INVALID|FE_DIVBYZERO),"undefined fit is not projected or normalized");
        if(symmetric_neighbors) {
            check(result.cell!=original&&result.relocated,"a symmetric connected structure uses the collision normal and topology");
            check(updater.diagnostics.unavoidable_splits==0,"undefined centroid does not force a neck break");
        } else {
            check(result.cell==original&&grid.material[original]==Gas&&grid.remaining(original)==0,"a true singleton is cleared");
            check(result.reaction_cells==0&&result.detached_si==4&&result.detached_cl==6,
                  "a singleton leaves as solid fragment inventory rather than chemical residual");
        }
        check(grid.silicon_inventory()+result.detached_si==si&&grid.chlorine_inventory()+result.detached_cl==cl,"Si and bound Cl close for connected and isolated updates");
        if(symmetric_neighbors)check(grid.remaining(original)>=26&&grid.remaining(original)<=30,"the surviving branch retains donor inventory plus conserved residual states");
    }
}


int main() {
    std::cout<<std::unitbuf;
    std::vector<std::pair<std::string,std::function<void()>>> tests={
        {"voxel storage and aperture",voxel_storage_and_hole},
        {"local quadratic surface geometry",quadratic_patch},
        {"normal projection",normal_projection},{"six coefficient least squares",six_coefficient_fit},
        {"fitted plane and y periodicity",fitted_plane_and_periodic_mapping},{"rational periodic orbit",closed_rational_orbit},{"y reflection and reverse path",reflection_in_y_plane},
        {"corner normal mapping",corner_normal_mapping},{"26-neighbor inventory and cache",neighbor_transfer_and_cache},
        {"solid-angle source distributions",source_solid_angle_distributions},{"shared engine inventory",shared_engine_inventory},
        {"actual collision projection and periodic image",collision_projection_and_periodic_image},
        {"multi-source connectivity equivalence",connectivity_matches_global_components},
        {"unbiased choices and capacity support",unbiased_candidate_ties_and_capacity_support},
        {"conservative deferred geometry",conservative_update_defers_geometry},
        {"fragment inventory and fixed anchors",fragment_release_and_anchors},
        {"collision normal at symmetric branch",collision_normal_at_symmetric_branch},
        {"connected inventory donor selection",connected_donor_selection},
        {"undefined normal removal",undefined_normal_removal},
        {"centroid degeneracy",centroid_degeneracy}
    };
    int failed=0;for(auto& [name,test]:tests){try{test();std::cout<<"PASS "<<name<<'\n';}
        catch(const std::exception& error){++failed;std::cerr<<"FAIL "<<name<<": "<<error.what()<<'\n';}}
    std::cout<<tests.size()-failed<<'/'<<tests.size()<<" spatial test groups passed\n";return failed?1:0;
}
