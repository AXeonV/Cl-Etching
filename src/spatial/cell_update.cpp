#include "cell_update.hpp"
#include <algorithm>
namespace etch::spatial {
CellRemoval CellUpdater::remove(int original) {
    Cell p=grid.coordinates(original);Patch patch=surface.patch(p.x,p.y,p.z);
    auto center=patch.project(patch.origin);Vec normal=patch.normal(center.u,center.v);
    double best=-center.distance;int selected=original;
    for(int dz=-1;dz<=1;++dz)for(int dy=-1;dy<=1;++dy)for(int dx=-1;dx<=1;++dx) {
        int z=p.z+dz;if(z<0||z>=grid.geometry.nz)continue;
        int cell=grid.index(p.x+dx,p.y+dy,z);
        if(cell==original||grid.material[cell]!=Silicon||!grid.boundary[cell])continue;
        Vec offset{double(dx),double(dy),double(dz)};
        double angle=std::acos(std::clamp(dot(unit(offset),normal),-1.0,1.0));
        double distance=patch.project(patch.origin+offset).distance;
        double score=normal_weight*angle-distance;
        if(score<best){best=score;selected=cell;}
    }
    Sites residual=grid.sites[original];
    if(selected!=original)grid.sites[original]=grid.sites[selected];
    grid.sites[selected]={};grid.remove_empty_cell(selected);
    Cell removed=grid.coordinates(selected);surface.invalidate(removed.x,removed.y,removed.z);
    return {residual,selected!=original};
}
}
