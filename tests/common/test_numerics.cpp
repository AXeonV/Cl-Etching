#include "common/implicit.hpp"
#include "common/periodic.hpp"
#include "common/traversal.hpp"
#include "common/vector.hpp"
#include <set>
#include <iostream>
#include <stdexcept>
#include <random>
using namespace etch;
void check(bool condition,const char* message){if(!condition)throw std::runtime_error(message);}
void close(double a,double b,double tolerance=1e-10){check(std::abs(a-b)<tolerance,"numerical mismatch");}
void bounded_roots() {
    std::vector<double> polynomial{1};
    for(double root:{0.1,0.3,0.6,0.9}){
        std::vector<double> next(polynomial.size()+1);
        for(size_t i=0;i<polynomial.size();++i){next[i]-=root*polynomial[i];next[i+1]+=polynomial[i];}polynomial=next;
    }
    auto endpoint=segment_roots({1+2e-16,-1});check(endpoint.size()==1,"retain a shared-face root within roundoff");close(endpoint[0],1);
    auto roots=segment_roots(polynomial);check(roots.size()==4,"all four crossings are isolated");
    double expected[]={0.1,0.3,0.6,0.9};for(int i=0;i<4;++i)close(roots[i],expected[i]);
    roots=segment_roots({0.25,-1,1});check(roots.size()==1,"retain a repeated tangent root");close(roots[0],0.5);
}
void quadratic_reproduction() {
    using Field=ImplicitCell<3>;std::array<Field::Node,8> nodes{};
    for(int corner=0;corner<8;++corner){double x=corner&1,y=(corner>>1)&1,z=(corner>>2)&1;
        nodes[corner]={z-0.2*x*x-0.1*y*y-0.6,-0.4*x,-0.2*y,0,1,0,0,0};}
    Field field(nodes,1);
    for(int i=0;i<100;++i){double x=(i+0.3)/100,y=(i+0.6)/100,z=(i+0.1)/100;
        close(field.value({x,y,z}),z-0.2*x*x-0.1*y*y-0.6);auto g=field.gradient({x,y,z});close(g[0],-0.4*x);close(g[1],-0.2*y);close(g[2],1);}
    auto roots=field.intersections({0.1,0.2,0},{0.8,0.7,1});check(roots.size()==1,"quadratic surface has one forward intersection");
    double t=roots[0];close(t-0.2*std::pow(0.1+0.7*t,2)-0.1*std::pow(0.2+0.5*t,2)-0.6,0);
}
ImplicitCell<3> adjacent_field(int offset) {
    using Field=ImplicitCell<3>;std::array<Field::Node,8> nodes{};
    for(int corner=0;corner<8;++corner){double x=offset+(corner&1),y=(corner>>1)&1,z=(corner>>2)&1;
        nodes[corner]={std::sin(x)+std::cos(y)+z*z+0.3*x*y+0.2*y*z+0.1*x*z,
                       std::cos(x)+0.3*y+0.1*z,-std::sin(y)+0.3*x+0.2*z,0.3,2*z+0.2*y+0.1*x,0.1,0.2,0};}
    return Field(nodes,1);
}
void face_continuity_and_bounds() {
    auto left=adjacent_field(0),right=adjacent_field(1);
    for(int i=0;i<100;++i){double y=(i+0.3)/100,z=(i+0.8)/100;
        close(left.value({1,y,z}),right.value({0,y,z}));auto a=left.gradient({1,y,z}),b=right.gradient({0,y,z});
        for(int axis=0;axis<3;++axis)close(a[axis],b[axis]);
        auto range=left.last_axis_range(z);double value=left.value({y,1-y,z});
        check(value>=range.first-1e-12&&value<=range.second+1e-12,"section Bernstein bounds contain values");
        check(value>=left.minimum-1e-12&&value<=left.maximum+1e-12,"volume Bernstein bounds contain values");}
}
void multiple_surface_crossings() {
    using Field=ImplicitCell<3>;std::array<Field::Node,8> nodes{};
    const double roots[3][3]={{0.1,0.4,0.8},{0.2,0.5,0.9},{0.3,0.6,0.7}};
    for(int corner=0;corner<8;++corner){double f[3],df[3];
        for(int axis=0;axis<3;++axis){double x=(corner>>axis)&1,a=x-roots[axis][0],b=x-roots[axis][1],c=x-roots[axis][2];f[axis]=a*b*c;df[axis]=a*b+a*c+b*c;}
        for(int mask=0;mask<8;++mask){double value=1;for(int axis=0;axis<3;++axis)value*=mask&(1<<axis)?df[axis]:f[axis];nodes[corner][mask]=value;}}
    Field field(nodes,1);auto found=field.intersections({0,0,0},{1,1,1});check(found.size()==9,"isolate all nine tricubic ray roots");
    for(int i=0;i<9;++i)close(found[i],0.1*(i+1),2e-7);
}
void periodic_closure() {
    close(horizontal_ray_period(1,2,32,32),32);
    close(horizontal_ray_period(2,3,30,28),420);
    close(horizontal_ray_period(0,-2,32,12),6);
    double dx=1/std::sqrt(5.0);close(horizontal_ray_period(dx,2*dx,32,32),32*std::sqrt(5.0));
}
void partition_derivatives() {
    ImplicitBlend<3> blend;ImplicitCell<3>::Node local{1.3,0.4,-0.2,0.3,0.8,0.1,0.2,0};
    blend.add({0.2,-0.3,0.4},{2,2,2},local);blend.add({-0.8,0.7,-0.6},{2,2,2},local);
    auto result=blend.normalized();for(int mask=0;mask<8;++mask)close(result[mask],local[mask]);
}
void vectors() {
    Vec3 x{1,0,0},y{0,1,0};close(norm(cross(x,y)-Vec3{0,0,1}),0);
    Vec3 d=unit({1,2,3}),n=unit({2,-1,1}),r=reflected(d,n);
    close(norm(r),1);close(dot(r,n),-dot(d,n));close(norm(reflected(r,n)-d),0);
}
void supercover_vertices_and_edges() {
    GridWalk<3> walk({0.5,0.5,0.5},unit({1,1,1}));std::set<std::array<int,3>> visited;
    for(int i=0;i<7;++i){walk.advance();visited.insert(walk.cell);}
    check(visited.size()==7 && visited.count({1,1,1}),"visit all seven vertex-touching cells");
    for(int mask=1;mask<8;++mask)check(visited.count({mask&1,(mask>>1)&1,(mask>>2)&1}),"vertex supercover subset");
    GridWalk<3> edge({0.5,0.5,0.5},unit({1,1,0}));visited.clear();
    for(int i=0;i<3;++i){edge.advance();visited.insert(edge.cell);}
    check(visited.size()==3 && visited.count({1,1,0}),"edge crossing visits three cells");
}
void planar_corner_supercover() {
    GridWalk<2> walk(std::array<double,2>{0.5,0.5},std::array<double,2>{1,1});
    std::array<std::array<int,2>,3> expected{{{1,0},{0,1},{1,1}}};
    for(int i=0;i<3;++i) {
        walk.advance();check(walk.cell==expected[i],"planar corner cells preserve traversal order");
        close(walk.entry,0.5);close(walk.exit(),i<2?0.5:1.5);
    }
    GridWalk<2> negative(std::array<double,2>{-0.5,0.5},std::array<double,2>{-1,0});
    negative.advance();check(negative.cell==std::array<int,2>{-2,0},"negative coordinates and zero direction component");
    close(negative.entry,0.5);close(negative.exit(),1.5);
}
void ray_box_clipping() {
    using Point=std::array<double,2>;
    auto plane=ray_box_interval<2>(Point{0,0},Point{0.5,0},Point{2,0},1);
    check(plane.has_value(),"parallel ray inside the box slab");close(plane->first,2);close(plane->last,6);
    check(!ray_box_interval<2>(Point{0,0},Point{0.5,0},Point{2,3},1),"parallel ray outside the slab");
    auto volume=ray_box_interval<3>(Vec3{4,0,0},Vec3{-2,0,0},Vec3{0,0,0},1);
    check(volume.has_value(),"negative ray direction intersects the box");close(volume->first,1.5);close(volume->last,2.5);
}
int main(){int failed=0;
    for(auto test:{bounded_roots,quadratic_reproduction,face_continuity_and_bounds,multiple_surface_crossings,periodic_closure,partition_derivatives,vectors,supercover_vertices_and_edges,planar_corner_supercover,ray_box_clipping}) {
        try{test();std::cout<<"PASS shared numerical check\n";}catch(const std::exception& error){++failed;std::cerr<<"FAIL "<<error.what()<<'\n';}}
    std::cout<<10-failed<<"/10 shared numerical groups passed\n";return failed?1:0;
}
