#include "cell_update.hpp"
#include <algorithm>

namespace etch::spatial {
CellUpdater::CellUpdater(Grid& g,Surface& s,double weight,unsigned long long seed)
    :grid(g),surface(s),normal_weight(weight),update(g,diagnostics,seed) {}

void CellUpdater::erase(int cell) {
    grid.remove_empty_cell(cell);
    Cell p=grid.coordinates(cell);surface.invalidate(p.x,p.y,p.z);
}
CellRemoval CellUpdater::remove(int original,Vec3 collision_normal,Vec3 collision_position) {
    auto erase_cell=[&](int cell){erase(cell);};
    if(update.isolated(original))return update.detach(original,erase_cell);
    Cell p=grid.coordinates(original);
    Patch patch=surface.patch(p.x,p.y,p.z);
    Vec3 normal=collision_normal;bool use_patch=patch.defined;double orientation=1;
    if(!patch.defined)++diagnostics.undefined_fit_removals;
    if(use_patch) {
        collision_position.x=patch.origin.x+std::remainder(collision_position.x-patch.origin.x,double(grid.geometry.nx));
        collision_position.y=patch.origin.y+std::remainder(collision_position.y-patch.origin.y,double(grid.geometry.ny));
        auto projection=patch.project(collision_position);
        auto fitted=patch.normal(projection.u,projection.v);double alignment=dot(fitted,collision_normal);
        if(alignment==0){use_patch=false;++diagnostics.tangent_fallbacks;}
        else {
            orientation=alignment>0?1:-1;diagnostics.reoriented_fits+=orientation<0;
            normal=fitted*orientation;
        }
    }
    double origin_distance=use_patch?orientation*patch.project(patch.origin).distance:0;
    constexpr int radius=1;
    std::vector<int> neighborhood;
    for(int dz=-radius;dz<=radius;++dz)for(int dy=-radius;dy<=radius;++dy)for(int dx=-radius;dx<=radius;++dx)
        if(p.z+dz>=0&&p.z+dz<grid.geometry.nz)neighborhood.push_back(grid.index(p.x+dx,p.y+dy,p.z+dz));
    auto branch=update.connectivity.local_silicon_branch(original,neighborhood);
    std::vector<RemovalCandidate> candidates{{original,-origin_distance,0}};
    for(int dz=-radius;dz<=radius;++dz)for(int dy=-radius;dy<=radius;++dy)for(int dx=-radius;dx<=radius;++dx) {
        int z=p.z+dz;if(z<0||z>=grid.geometry.nz)continue;
        int cell=grid.index(p.x+dx,p.y+dy,z);
        if(cell==original||grid.material[cell]!=Silicon||!grid.boundary[cell])continue;
        if(!branch.contains(cell)){++diagnostics.disconnected_donors;continue;}
        Vec3 offset{double(dx),double(dy),double(dz)};
        double angle=std::acos(std::clamp(dot(unit(offset),normal),-1.0,1.0));
        double distance=use_patch?orientation*patch.project(patch.origin+offset).distance:dot(offset,collision_normal);
        candidates.push_back({cell,normal_weight*angle-distance,dot(offset,offset)});
    }
    return update.apply(original,{std::move(branch),std::move(candidates)},erase_cell);
}
}
