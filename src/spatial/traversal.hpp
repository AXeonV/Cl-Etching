#pragma once
#include "vector.hpp"
#include <array>
#include <algorithm>
#include <limits>
namespace etch::spatial {
// Supercover DDA: face crossings visit one cell, edges three, vertices seven.
class VoxelWalk {
    std::array<int,3> step;
    std::array<double,3> next,delta;
    std::array<std::array<int,3>,7> touched;
    int touched_count=0,touched_index=0;
public:
    std::array<int,3> cell;
    double entry=0;
    VoxelWalk(Vec p,Vec d) {
        for(int axis=0;axis<3;++axis) {
            cell[axis]=static_cast<int>(std::floor(p[axis]));step[axis]=d[axis]>0?1:-1;
            next[axis]=d[axis]==0?std::numeric_limits<double>::infinity():
                (cell[axis]+(step[axis]>0)-p[axis])/d[axis];
            delta[axis]=d[axis]==0?std::numeric_limits<double>::infinity():1/std::abs(d[axis]);
        }
    }
    double exit() const {
        return touched_index+1<touched_count?entry:std::min({next[0],next[1],next[2]});
    }
    void advance() {
        if(touched_index+1<touched_count) { cell=touched[++touched_index];return; }
        entry=std::min({next[0],next[1],next[2]});int mask=0;
        for(int axis=0;axis<3;++axis)if(next[axis]==entry)mask|=1<<axis;
        touched_count=0;
        for(int subset=1;subset<=mask;++subset)if((subset&mask)==subset) {
            auto candidate=cell;
            for(int axis=0;axis<3;++axis)if(subset&(1<<axis))candidate[axis]+=step[axis];
            touched[touched_count++]=candidate;
        }
        for(int axis=0;axis<3;++axis)if(mask&(1<<axis))next[axis]+=delta[axis];
        touched_index=0;cell=touched[0];
    }
};
}
