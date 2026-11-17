#include "grid.hpp"
#include <algorithm>
namespace etch::spatial {
Grid::Grid(Geometry g):geometry(g),material(g.nx*g.ny*g.nz),boundary(material.size()),sites(material.size()) {
    int left=(g.nx-g.opening)/2,front=(g.ny-g.opening_y)/2;
    for(int z=0;z<g.nz;++z)for(int y=0;y<g.ny;++y)for(int x=0;x<g.nx;++x) {
        int cell=index(x,y,z);
        if(z>=g.surface_z) { material[cell]=Silicon;sites[cell][0]=g.atoms_per_cell; }
        else if(z>=g.surface_z-g.mask_thickness &&
                (x<left||x>=left+g.opening||y<front||y>=front+g.opening_y)) material[cell]=Mask;
    }
    refresh(g.nx/2,g.ny/2,g.nz/2,std::max({g.nx,g.ny,g.nz}));
}
std::array<int,6> Grid::face_neighbors(int cell) const {
    Cell p=coordinates(cell);
    return {index(p.x-1,p.y,p.z),index(p.x+1,p.y,p.z),index(p.x,p.y-1,p.z),index(p.x,p.y+1,p.z),
            p.z>0?index(p.x,p.y,p.z-1):-1,p.z+1<geometry.nz?index(p.x,p.y,p.z+1):-1};
}
Material Grid::at(int x,int y,int z) const {
    if(z<0)return Gas;
    if(z>=geometry.nz)return Silicon;
    return static_cast<Material>(material[index(x,y,z)]);
}
void Grid::wrap_position(Vec3& p) const {
    p.x-=geometry.nx*std::floor(p.x/geometry.nx);
    p.y-=geometry.ny*std::floor(p.y/geometry.ny);
}
void Grid::refresh(int x,int y,int z,int radius) {
    int xmin=x-std::min(radius,geometry.nx/2),ymin=y-std::min(radius,geometry.ny/2);
    int xmax=radius>=geometry.nx/2?xmin+geometry.nx-1:x+radius;
    int ymax=radius>=geometry.ny/2?ymin+geometry.ny-1:y+radius;
    for(int zz=std::max(0,z-radius);zz<=std::min(geometry.nz-1,z+radius);++zz)
        for(int yy=ymin;yy<=ymax;++yy)for(int xx=xmin;xx<=xmax;++xx) {
            int cell=index(xx,yy,zz);
            boundary[cell]=at(xx,yy,zz)!=Gas &&
                (at(xx-1,yy,zz)==Gas||at(xx+1,yy,zz)==Gas||at(xx,yy-1,zz)==Gas||
                 at(xx,yy+1,zz)==Gas||at(xx,yy,zz-1)==Gas||at(xx,yy,zz+1)==Gas);
            if(boundary[cell])surface_cells.insert(cell);else surface_cells.erase(cell);
        }
}
void Grid::remove_empty_cell(int cell) {
    material[cell]=Gas;Cell p=coordinates(cell);refresh(p.x,p.y,p.z,2);
}
}
