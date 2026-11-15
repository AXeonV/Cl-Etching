#include "common/implicit.hpp"
#include "common/periodic.hpp"
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
int main(){int failed=0;
    for(auto test:{bounded_roots,quadratic_reproduction,face_continuity_and_bounds,multiple_surface_crossings,periodic_closure,partition_derivatives}) {
        try{test();std::cout<<"PASS shared numerical check\n";}catch(const std::exception& error){++failed;std::cerr<<"FAIL "<<error.what()<<'\n';}}
    std::cout<<6-failed<<"/6 shared numerical groups passed\n";return failed?1:0;
}
