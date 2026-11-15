#pragma once
#include <cmath>
#include <algorithm>
#include <cstdint>
#include <numeric>
#include <limits>
namespace etch {
// A horizontal ray in a rectangular torus closes when both lateral displacements
// are integer periods. Reduce the exact binary significands; no angular tolerance.
inline double horizontal_ray_period(double dx,double dy,int nx,int ny) {
    if(dx==0)return ny/std::abs(dy);
    if(dy==0)return nx/std::abs(dx);
    auto parts=[](double value,int& exponent) {
        double mantissa=std::frexp(std::abs(value),&exponent);
        uint64_t significant=static_cast<uint64_t>(std::ldexp(mantissa,53));exponent-=53;
        while((significant&1)==0){significant>>=1;++exponent;}return significant;
    };
    int ex,ey;uint64_t sx=parts(dx,ex),sy=parts(dy,ey),x=nx,y=ny;
    auto cancel=[](uint64_t& a,uint64_t& b){uint64_t g=std::gcd(a,b);a/=g;b/=g;};
    cancel(sx,sy);cancel(x,y);cancel(sx,x);cancel(sy,y);
    int shift=ex-ey;
    while(shift>0&&(x&1)==0){x>>=1;--shift;}
    while(shift<0&&(y&1)==0){y>>=1;++shift;}
    long double cycles=std::ldexp(static_cast<long double>(sx)*y,std::max(shift,0));
    return static_cast<double>(cycles*nx/std::abs(dx));
}
}
