#include "backend.hpp"
namespace etch::spatial {
void Backend::write_surface(Surface& surface,const Grid& grid,const std::string& path) {
    std::ofstream file(path);file<<"x,y,z,ux,uy,uz,vx,vy,vz,nx,ny,nz,a,b,c,d,e,f,cell\n"<<std::setprecision(15);
    for(int cell:grid.surface_cells) {
        Cell p=grid.coordinates(cell);const auto& s=surface.patch(p.x,p.y,p.z);
        if(!s.defined)continue;
        file<<s.origin.x<<','<<s.origin.y<<','<<s.origin.z<<','<<s.u_axis.x<<','<<s.u_axis.y<<','<<s.u_axis.z<<','
            <<s.v_axis.x<<','<<s.v_axis.y<<','<<s.v_axis.z<<','<<s.outward.x<<','<<s.outward.y<<','<<s.outward.z<<','
            <<s.a<<','<<s.b<<','<<s.c<<','<<s.d<<','<<s.e<<','<<s.f<<','<<cell<<'\n';
    }
}
void Backend::write_field(Surface& surface,const Grid& grid,const std::string& path) {
    std::ofstream file(path,std::ios::binary);int n=surface.resolution();
    for(int z=0;z<=grid.geometry.nz*n;++z)for(int y=0;y<=grid.geometry.ny*n;++y)for(int x=0;x<=grid.geometry.nx*n;++x) {
        auto node=surface.node(x,y,z);file.write(reinterpret_cast<const char*>(node.data()),sizeof(node));
    }
}
void Backend::write_geometry(Summary& out,const Geometry& g) {
    out.field("dimension",3);out.field("nx",g.nx);out.field("ny",g.ny);out.field("nz",g.nz);
    out.field("cell_nm",g.cell_nm);out.field("surface_z",g.surface_z);
    out.field("opening_cells",g.opening);out.field("opening_y_cells",g.opening_y);out.field("mask_cells",g.mask_thickness);
}
}
