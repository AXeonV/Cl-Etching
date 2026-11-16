#pragma once
#include <cmath>
#include <optional>
#include <array>
#include <algorithm>
#include <limits>
namespace etch {
// Visit every cell touched at a grid crossing, including corner and edge neighbors.
template<int D> class GridWalk {
    std::array<int,D> step;
    std::array<double,D> next,delta;
    std::array<std::array<int,D>,(1<<D)-1> touched;
    int touched_count=0,touched_index=0;
public:
    std::array<int,D> cell;
    double entry=0;
    template<class Vector> GridWalk(Vector p,Vector d) {
        for(int axis=0;axis<D;++axis) {
            cell[axis]=static_cast<int>(std::floor(p[axis]));step[axis]=d[axis]>0?1:-1;
            next[axis]=d[axis]==0?std::numeric_limits<double>::infinity():
                (cell[axis]+(step[axis]>0)-p[axis])/d[axis];
            delta[axis]=d[axis]==0?std::numeric_limits<double>::infinity():1/std::abs(d[axis]);
        }
    }
    double exit() const {
        return touched_index+1<touched_count?entry:*std::min_element(next.begin(),next.end());
    }
    void advance() {
        if(touched_index+1<touched_count) { cell=touched[++touched_index];return; }
        entry=*std::min_element(next.begin(),next.end());int mask=0;
        for(int axis=0;axis<D;++axis)if(next[axis]==entry)mask|=1<<axis;
        touched_count=0;
        for(int subset=1;subset<=mask;++subset)if((subset&mask)==subset) {
            auto candidate=cell;
            for(int axis=0;axis<D;++axis)if(subset&(1<<axis))candidate[axis]+=step[axis];
            touched[touched_count++]=candidate;
        }
        for(int axis=0;axis<D;++axis)if(mask&(1<<axis))next[axis]+=delta[axis];
        touched_index=0;cell=touched[0];
    }
};
struct RayInterval { double first=0,last=std::numeric_limits<double>::infinity(); };
template<int D,class Vector>
std::optional<RayInterval> ray_box_interval(Vector point,Vector direction,Vector center,double half_width) {
    RayInterval interval;
    for(int axis=0;axis<D;++axis) {
        if(direction[axis]==0) {
            if(point[axis]<center[axis]-half_width||point[axis]>center[axis]+half_width)return std::nullopt;
        } else {
            double first=(center[axis]-half_width-point[axis])/direction[axis];
            double last=(center[axis]+half_width-point[axis])/direction[axis];
            if(first>last)std::swap(first,last);
            interval.first=std::max(interval.first,first);interval.last=std::min(interval.last,last);
        }
    }
    if(interval.last<interval.first)return std::nullopt;
    return interval;
}
}
