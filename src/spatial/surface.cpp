#include "surface.hpp"
#include "common/linear.hpp"
#include <algorithm>
namespace etch::spatial {
Vec Patch::normal(double u,double v) const {
    return unit(outward-u_axis*(2*a*u+c*v+d)-v_axis*(2*b*v+c*u+e));
}
Projection Patch::project(Vec point) const {
    Vec q=point-origin;double u0=dot(q,u_axis),v0=dot(q,v_axis),w0=dot(q,outward);
    double bound=std::abs(height(u0,v0)-w0),best_u=u0,best_v=v0,best=bound*bound;
    auto objective=[&](double u,double v) {
        double h=height(u,v)-w0;
        return (u-u0)*(u-u0)+(v-v0)*(v-v0)+h*h;
    };
    // Multi-start Newton normal projection, with a descent line search.
    for(int su=-1;su<=1;++su)for(int sv=-1;sv<=1;++sv) {
        double u=u0+su*bound,v=v0+sv*bound;
        for(int iteration=0;iteration<40;++iteration) {
            double h=height(u,v)-w0,fu=2*a*u+c*v+d,fv=2*b*v+c*u+e;
            double gu=u-u0+h*fu,gv=v-v0+h*fv;
            if(gu*gu+gv*gv<1e-20)break;
            double h11=1+fu*fu+2*a*h,h22=1+fv*fv+2*b*h,h12=fu*fv+c*h;
            double determinant=h11*h22-h12*h12,du=gu,dv=gv;
            if(h11>0 && determinant>0) {
                du=(h22*gu-h12*gv)/determinant;dv=(h11*gv-h12*gu)/determinant;
            }
            double value=objective(u,v),scale=1;bool advanced=false;
            for(int step=0;step<24;++step) {
                double nu=u-scale*du,nv=v-scale*dv;
                if(objective(nu,nv)<value) { u=nu;v=nv;advanced=true;break; }
                scale*=0.5;
            }
            if(!advanced)break;
        }
        double value=objective(u,v);
        if(value<best) { best=value;best_u=u;best_v=v; }
    }
    Vec delta=point-position(best_u,best_v);
    return {best_u,best_v,std::copysign(std::sqrt(best),dot(delta,normal(best_u,best_v)))};
}
Surface::Surface(Grid& g,int r,int n):grid(g),radius(r),subdivisions(n),nearby_cells(g.material.size()),dirty_candidates(g.material.size(),1) {
    node_values.resize(g.geometry.nx*n*g.geometry.ny*n*(g.geometry.nz*n+1));dirty_nodes.assign(node_values.size(),1);
}
const Patch& Surface::patch(int x,int y,int z) {
    int cell=grid.index(x,y,z);auto found=patches.find(cell);
    if(found!=patches.end())return found->second;
    std::vector<Vec> points;Vec centroid{};int count=0;
    for(int dz=-radius;dz<=radius;++dz)for(int dy=-radius;dy<=radius;++dy)for(int dx=-radius;dx<=radius;++dx) {
        int zz=z+dz;if(zz<0||zz>=grid.geometry.nz)continue;
        int id=grid.index(x+dx,y+dy,zz);Vec offset{double(dx),double(dy),double(dz)};
        if(grid.material[id]!=Gas){centroid=centroid+offset;++count;}
        if(grid.boundary[id])points.push_back(offset);
    }
    centroid=centroid*(1.0/count);
    Patch result;result.origin={grid.wrap_x(x)+0.5,grid.wrap_y(y)+0.5,z+0.5};
    ++fit_count;
    if(norm(centroid)==0) { result.defined=false;++undefined_normal_fits;return patches.emplace(cell,result).first->second; }
    result.outward=unit(centroid*-1);
    Vec n=result.outward,reference;
    if(std::abs(n.x)<=std::abs(n.y)&&std::abs(n.x)<=std::abs(n.z))reference={1,0,0};
    else if(std::abs(n.y)<=std::abs(n.z))reference={0,1,0};else reference={0,0,1};
    result.u_axis=unit(cross(reference,n));result.v_axis=cross(n,result.u_axis);
    double matrix[6][6]={},rhs[6]={},coeff[6]={};
    for(Vec p:points) {
        double u=dot(p,result.u_axis),v=dot(p,result.v_axis),w=dot(p,n);
        double row[]={u*u,v*v,u*v,u,v,1};
        for(int i=0;i<6;++i){rhs[i]+=row[i]*w;for(int j=0;j<6;++j)matrix[i][j]+=row[i]*row[j];}
    }
    if(least_squares(matrix,rhs,coeff)<6)++rank_deficient_fits;
    result.a=coeff[0];result.b=coeff[1];result.c=coeff[2];
    result.d=coeff[3];result.e=coeff[4];result.f=coeff[5];
    return patches.emplace(cell,result).first->second;
}
const std::vector<int>& Surface::candidates(int x,int y,int z) {
    int id=grid.index(x,y,z);if(!dirty_candidates[id])return nearby_cells[id];
    auto& cells=nearby_cells[id];cells.clear();int reach=radius+1;
    int xmin=x-std::min(reach,grid.geometry.nx/2),ymin=y-std::min(reach,grid.geometry.ny/2);
    int xmax=reach>=grid.geometry.nx/2?xmin+grid.geometry.nx-1:x+reach;
    int ymax=reach>=grid.geometry.ny/2?ymin+grid.geometry.ny-1:y+reach;
    for(int zz=std::max(0,z-reach);zz<=std::min(grid.geometry.nz-1,z+reach);++zz)
        for(int yy=ymin;yy<=ymax;++yy)for(int xx=xmin;xx<=xmax;++xx)
            if(grid.boundary[grid.index(xx,yy,zz)])cells.push_back(grid.index(xx,yy,zz));
    dirty_candidates[id]=0;return cells;
}
bool Surface::can_hit_horizontal_ray(double z) {
    double level=z*subdivisions;int cell=static_cast<int>(std::floor(level));
    if(cell<0||cell>=grid.geometry.nz*subdivisions)return false;
    for(int y=0;y<grid.geometry.ny*subdivisions;++y)for(int x=0;x<grid.geometry.nx*subdivisions;++x) {
        auto [low,high]=field(x,y,cell).last_axis_range(level-cell);
        if(low<0&&high>0)return true;
    }
    return false;
}

Surface::Field::Node Surface::node(int x,int y,int z) {
    int width=grid.geometry.nx*subdivisions,height=grid.geometry.ny*subdivisions;
    int id=(z*height+(y%height+height)%height)*width+(x%width+width)%width;
    if(!dirty_nodes[id])return node_values[id];
    Vec point{double(x)/subdivisions,double(y)/subdivisions,double(z)/subdivisions};ImplicitBlend<3> blend;
    double support=radius+0.5;
    std::array<double,3> window{std::min(support,grid.geometry.nx/2.0),std::min(support,grid.geometry.ny/2.0),support};
    for(int cell:candidates(int(std::floor(point.x)),int(std::floor(point.y)),std::min(int(std::floor(point.z)),grid.geometry.nz-1))) {
        Cell p=grid.coordinates(cell);Patch fit=patch(p.x,p.y,p.z);if(!fit.defined)continue;
        fit.origin.x+=grid.geometry.nx*std::round((point.x-fit.origin.x)/grid.geometry.nx);
        fit.origin.y+=grid.geometry.ny*std::round((point.y-fit.origin.y)/grid.geometry.ny);
        Vec relative=point-fit.origin;double u=dot(relative,fit.u_axis),v=dot(relative,fit.v_axis),w=dot(relative,fit.outward);
        double fu=2*fit.a*u+fit.c*v+fit.d,fv=2*fit.b*v+fit.c*u+fit.e;
        Vec gradient=fit.outward-fit.u_axis*fu-fit.v_axis*fv;
        auto mixed=[&](int i,int j){return -2*fit.a*fit.u_axis[i]*fit.u_axis[j]-2*fit.b*fit.v_axis[i]*fit.v_axis[j]-fit.c*(fit.u_axis[i]*fit.v_axis[j]+fit.v_axis[i]*fit.u_axis[j]);};
        Field::Node local{w-fit.height(u,v),gradient.x,gradient.y,mixed(0,1),gradient.z,mixed(0,2),mixed(1,2),0};
        blend.add({relative.x,relative.y,relative.z},window,local);
    }
    Field::Node value{};
    if(blend.weights[0]>0)value=blend.normalized();
    else value[0]=grid.at(int(std::floor(point.x)),int(std::floor(point.y)),int(std::floor(point.z)))==Gas?radius+1:-(radius+1);
    dirty_nodes[id]=0;return node_values[id]=value;
}

const Surface::Field& Surface::field(int x,int y,int z) {
    int width=grid.geometry.nx*subdivisions,height=grid.geometry.ny*subdivisions;
    int id=(z*height+(y%height+height)%height)*width+(x%width+width)%width;
    auto found=fields.find(id);if(found!=fields.end())return found->second;
    std::array<Field::Node,8> samples;
    for(int corner=0;corner<8;++corner)samples[corner]=node(x+(corner&1),y+((corner>>1)&1),z+((corner>>2)&1));
    return fields.emplace(id,Field(samples,1.0/subdivisions)).first->second;
}
void Surface::invalidate(int x,int y,int z) {
    int reach=(2*radius+3)*subdivisions+1,width=grid.geometry.nx*subdivisions,height=grid.geometry.ny*subdivisions;
    for(int dz=-reach;dz<=reach;++dz)for(int dy=-reach;dy<=reach;++dy)for(int dx=-reach;dx<=reach;++dx) {
        int zz=z*subdivisions+dz,xx=x*subdivisions+dx,yy=y*subdivisions+dy;
        if(zz>=0&&zz<=grid.geometry.nz*subdivisions){int id=(zz*height+(yy%height+height)%height)*width+(xx%width+width)%width;dirty_nodes[id]=1;fields.erase(id);}
    }
    for(int dz=-radius-2;dz<=radius+2;++dz)for(int dy=-radius-2;dy<=radius+2;++dy)for(int dx=-radius-2;dx<=radius+2;++dx) {
        int zz=z+dz;if(zz<0||zz>=grid.geometry.nz)continue;
        int cell=grid.index(x+dx,y+dy,zz);dirty_candidates[cell]=1;
        if(std::abs(dx)<=radius+1&&std::abs(dy)<=radius+1&&std::abs(dz)<=radius+1)patches.erase(cell);
    }
}
}
