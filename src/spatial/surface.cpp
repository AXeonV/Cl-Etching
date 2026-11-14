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
std::optional<double> Patch::intersection(Vec p,Vec direction,double start,double end) const {
    Vec q=p-origin;
    double u=dot(q,u_axis),v=dot(q,v_axis),w=dot(q,outward);
    double du=dot(direction,u_axis),dv=dot(direction,v_axis),dw=dot(direction,outward);
    double A=a*du*du+b*dv*dv+c*du*dv;
    double B=2*a*u*du+2*b*v*dv+c*(u*dv+v*du)+d*du+e*dv-dw;
    double C=height(u,v)-w;
    double roots[2];int count=0;
    if(std::abs(A)<1e-12) { if(std::abs(B)>1e-12)roots[count++]=-C/B; }
    else {
        double discriminant=B*B-4*A*C;
        if(discriminant>=0) {
            double qroot=-0.5*(B+std::copysign(std::sqrt(discriminant),B));
            roots[count++]=qroot/A;
            if(std::abs(qroot)>1e-15)roots[count++]=C/qroot;
        }
    }
    double nearest=std::numeric_limits<double>::infinity();
    for(int i=0;i<count;++i) {
        double t=roots[i],uu=u+du*t,vv=v+dv*t;
        if(t>=start && t<=end && uu>=u_min && uu<=u_max && vv>=v_min && vv<=v_max
           && dot(direction,normal(uu,vv))<-1e-9)nearest=std::min(nearest,t);
    }
    if(std::isfinite(nearest))return nearest;
    return std::nullopt;
}
std::pair<double,double> Patch::z_bounds() const {
    double A=outward.z*a,B=outward.z*b,C=outward.z*c;
    double D=u_axis.z+outward.z*d,E=v_axis.z+outward.z*e,F=origin.z+outward.z*f;
    double low=std::numeric_limits<double>::infinity(),high=-low;
    auto sample=[&](double u,double v) {
        if(u<u_min||u>u_max||v<v_min||v>v_max)return;
        double value=A*u*u+B*v*v+C*u*v+D*u+E*v+F;
        low=std::min(low,value);high=std::max(high,value);
    };
    for(double u:{u_min,u_max})for(double v:{v_min,v_max})sample(u,v);
    if(B!=0)for(double u:{u_min,u_max})sample(u,-(C*u+E)/(2*B));
    if(A!=0)for(double v:{v_min,v_max})sample(-(C*v+D)/(2*A),v);
    double determinant=4*A*B-C*C;
    if(determinant!=0)sample((C*E-2*B*D)/determinant,(C*D-2*A*E)/determinant);
    return {low,high};
}

Surface::Surface(Grid& g,int r):grid(g),radius(r),nearby_cells(g.material.size()),dirty_candidates(g.material.size(),1) {}
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
    result.u_min=result.v_min=std::numeric_limits<double>::infinity();
    result.u_max=result.v_max=-std::numeric_limits<double>::infinity();
    for(Vec p:points) {
        double u=dot(p,result.u_axis),v=dot(p,result.v_axis),w=dot(p,n);
        double row[]={u*u,v*v,u*v,u,v,1};
        result.u_min=std::min(result.u_min,u);result.u_max=std::max(result.u_max,u);
        result.v_min=std::min(result.v_min,v);result.v_max=std::max(result.v_max,v);
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
bool Surface::can_hit_horizontal_ray(double z,Vec direction) {
    for(int cell:grid.surface_cells) {
        Cell p=grid.coordinates(cell);const Patch& s=patch(p.x,p.y,p.z);
        if(!s.defined)continue;
        auto [low,high]=s.z_bounds();
        if(z<low||z>high||low==high)continue;
        double minimum=std::numeric_limits<double>::infinity();
        for(double u:{s.u_min,s.u_max})for(double v:{s.v_min,s.v_max}) {
            Vec normal=s.outward-s.u_axis*(2*s.a*u+s.c*v+s.d)-s.v_axis*(2*s.b*v+s.c*u+s.e);
            minimum=std::min(minimum,dot(direction,normal));
        }
        if(minimum<-1e-9)return true;
    }
    return false;
}

void Surface::invalidate(int x,int y,int z) {
    for(int dz=-radius-2;dz<=radius+2;++dz)for(int dy=-radius-2;dy<=radius+2;++dy)for(int dx=-radius-2;dx<=radius+2;++dx) {
        int zz=z+dz;if(zz<0||zz>=grid.geometry.nz)continue;
        int cell=grid.index(x+dx,y+dy,zz);dirty_candidates[cell]=1;
        if(std::abs(dx)<=radius+1&&std::abs(dy)<=radius+1&&std::abs(dz)<=radius+1)patches.erase(cell);
    }
}
}
