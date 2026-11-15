#include "surface.hpp"
#include "common/linear.hpp"
#include <algorithm>
#include <cmath>

namespace etch::planar {
Vec Curve::position(double u) const {
    return origin + tangent * u + outward * (a * u * u + b * u + c);
}
Vec Curve::normal(double u) const { return unit(outward - tangent * (2 * a * u + b)); }
// The nearest point solves (u-u0) + (f(u)-v0)*f'(u) = 0.
CurveProjection Curve::project(Vec point) const {
    Vec local = point - origin;
    double u0 = dot(local, tangent), v0 = dot(local, outward);
    double u = (u0 + b * (v0 - c)) / (1 + b*b);
    if (a != 0) {
        auto slope_equation = [&](double x) {
            return (x-u0) + (a*x*x+b*x+c-v0)*(2*a*x+b);
        };
        // A closer point cannot lie farther horizontally than this known distance.
        double bound = std::abs(a*u0*u0+b*u0+c-v0);
        std::vector<double> intervals{u0-bound, u0+bound};
        double A = 6*a*a, B = 6*a*b, C = 1+b*b+2*a*(c-v0);
        double discriminant = B*B - 4*A*C;
        if (discriminant > 0) {
            double q = -0.5*(B + std::copysign(std::sqrt(discriminant), B));
            for (double critical : {q/A, C/q})
                if (critical > intervals[0] && critical < intervals[1]) intervals.push_back(critical);
        }
        std::sort(intervals.begin(), intervals.end());
        u = u0;
        double nearest = norm(point - position(u));
        for (size_t i = 1; i < intervals.size(); ++i) {
            double left = intervals[i-1], right = intervals[i];
            double gl = slope_equation(left), gr = slope_equation(right);
            if ((gl > 0 && gr > 0) || (gl < 0 && gr < 0)) continue;
            if (gl == 0) right = left;
            else if (gr == 0) left = right;
            for (int step = 0; step < 60; ++step) {
                double middle = (left+right)/2, gm = slope_equation(middle);
                if (gm == 0) { left = right = middle; break; }
                if ((gl <= 0 && gm <= 0) || (gl >= 0 && gm >= 0)) { left = middle; gl = gm; }
                else right = middle;
            }
            double candidate = (left+right)/2, distance = norm(point-position(candidate));
            if (distance < nearest) { nearest = distance; u = candidate; }
        }
    }
    Vec delta = point-position(u);
    return {u, std::copysign(norm(delta), dot(delta, normal(u)))};
}

Surface::Surface(Grid& g, int r,int n) : grid(g), radius(r),subdivisions(n), curves(g.material.size()), dirty(g.material.size(), 1),
    dirty_candidates(g.material.size(),1), nearby_cells(g.material.size()) {
    node_values.resize(g.geometry.nx*n*(g.geometry.nz*n+1));dirty_nodes.assign(node_values.size(),1);
}
const Curve& Surface::curve(int x, int z) {
    int id = grid.index(x, z);
    if (!dirty[id]) return curves[id];
    std::vector<Vec> points;
    Vec solid_centroid{};
    int solid_count = 0;
    for (int dz = -radius; dz <= radius; ++dz) for (int dx = -radius; dx <= radius; ++dx) {
        int zz = z+dz;
        if (zz < 0 || zz >= grid.geometry.nz) continue;
        int cell = grid.index(x+dx,zz);
        if (grid.material[cell] != Gas) {
            solid_centroid = solid_centroid+Vec{static_cast<double>(dx),static_cast<double>(dz)};
            ++solid_count;
        }
        if (grid.boundary[cell]) points.push_back({static_cast<double>(dx),static_cast<double>(dz)});
    }
    solid_centroid = solid_centroid*(1.0/solid_count);
    if (solid_centroid.x==0 && solid_centroid.z==0) {
        curves[id] = {}; curves[id].origin = {x+0.5,z+0.5}; curves[id].defined = false;
        dirty[id] = 0; ++fit_count; ++undefined_normal_fits;
        return curves[id];
    }
    Vec outward = unit(solid_centroid*-1);
    Vec tangent{-outward.z,outward.x};
    double matrix[3][3] = {}, rhs[3] = {}, coefficients[3] = {};
    for (Vec p : points) {
        double u = dot(p, tangent), v = dot(p, outward), row[] = {u*u, u, 1};
        for (int j = 0; j < 3; ++j) {
            rhs[j] += row[j] * v;
            for (int k = 0; k < 3; ++k) matrix[j][k] += row[j] * row[k];
        }
    }
    int rank = least_squares(matrix, rhs, coefficients);
    if (rank < 3) ++rank_deficient_fits;
    curves[id] = {{x+0.5,z+0.5},tangent,outward,
                  coefficients[0],coefficients[1],coefficients[2]};
    dirty[id] = 0;
    ++fit_count;
    return curves[id];
}
const std::vector<int>& Surface::candidates(int x, int z) {
    int id = grid.index(x,z);
    if (!dirty_candidates[id]) return nearby_cells[id];
    auto& cells = nearby_cells[id]; cells.clear();
    int reach = radius+1;
    for (int dz=-reach;dz<=reach;++dz) for(int dx=-reach;dx<=reach;++dx) {
        int zz=z+dz;
        if (zz>=0 && zz<grid.geometry.nz && grid.boundary[grid.index(x+dx,zz)])
            cells.push_back(grid.index(x+dx,zz));
    }
    dirty_candidates[id]=0;
    return cells;
}


Surface::Field::Node Surface::node(int x,int z) {
    int width=grid.geometry.nx*subdivisions,id=z*width+(x%width+width)%width;
    if(!dirty_nodes[id])return node_values[id];
    Vec point{double(x)/subdivisions,double(z)/subdivisions};ImplicitBlend<2> blend;
    double support=radius+0.5;std::array<double,2> window{std::min(support,grid.geometry.nx/2.0),support};
    for(int cell:candidates(int(std::floor(point.x)),std::min(int(std::floor(point.z)),grid.geometry.nz-1))) {
        Curve fit=curve(cell%grid.geometry.nx,cell/grid.geometry.nx);if(!fit.defined)continue;
        fit.origin.x+=grid.geometry.nx*std::round((point.x-fit.origin.x)/grid.geometry.nx);
        Vec relative=point-fit.origin;double u=dot(relative,fit.tangent),v=dot(relative,fit.outward),slope=2*fit.a*u+fit.b;
        Vec gradient=fit.outward-fit.tangent*slope;
        Field::Node local{v-(fit.a*u*u+fit.b*u+fit.c),gradient.x,gradient.z,-2*fit.a*fit.tangent.x*fit.tangent.z};
        blend.add({relative.x,relative.z},window,local);
    }
    Field::Node value{};
    if(blend.weights[0]>0)value=blend.normalized();
    else value[0]=grid.at(int(std::floor(point.x)),int(std::floor(point.z)))==Gas?radius+1:-(radius+1);
    dirty_nodes[id]=0;return node_values[id]=value;
}

const Surface::Field& Surface::field(int x,int z) {
    int width=grid.geometry.nx*subdivisions,id=z*width+(x%width+width)%width;
    auto found=fields.find(id);if(found!=fields.end())return found->second;
    std::array<Field::Node,4> samples;
    for(int corner=0;corner<4;++corner)samples[corner]=node(x+(corner&1),z+((corner>>1)&1));
    return fields.emplace(id,Field(samples,1.0/subdivisions)).first->second;
}
void Surface::invalidate(int x, int z) {
    int reach=(2*radius+3)*subdivisions+1,width=grid.geometry.nx*subdivisions;
    for(int dz=-reach;dz<=reach;++dz)for(int dx=-reach;dx<=reach;++dx) {
        int zz=z*subdivisions+dz,xx=x*subdivisions+dx;
        if(zz>=0&&zz<=grid.geometry.nz*subdivisions){int id=zz*width+(xx%width+width)%width;dirty_nodes[id]=1;fields.erase(id);}
    }
    for (int dz=-radius-2;dz<=radius+2;++dz) for(int dx=-radius-2;dx<=radius+2;++dx) {
        if (z+dz<0 || z+dz>=grid.geometry.nz) continue;
        int id=grid.index(x+dx,z+dz);
        dirty_candidates[id]=1;
        if (std::abs(dx)<=radius+1 && std::abs(dz)<=radius+1) dirty[id]=1;
    }
}
}
