#pragma once
#include "material.hpp"
#include "random.hpp"
#include <optional>

namespace etch {
// Copy the donor and redistribute the residual over equally probable free slots.
inline std::optional<CellRemoval> transfer_inventory(std::vector<Sites>& sites,int original,int selected,
        const std::vector<int>& recipients,int capacity,Random& rng) {
    Sites pending=sites[original];
    int residual=site_count(pending),available=0;
    for(int cell:recipients) {
        int source=cell==original?selected:cell;
        available+=capacity-site_count(sites[source]);
    }
    if(available<residual)return std::nullopt;

    if(selected!=original)sites[original]=sites[selected];
    sites[selected]={};
    CellRemoval removal{selected,selected!=original};
    removal.redistributed_si=residual;
    for(int state=1;state<4;++state)removal.redistributed_cl+=state*pending[state];

    std::vector<int> free_slots;
    for(int cell:recipients)free_slots.push_back(capacity-site_count(sites[cell]));
    while(residual>0) {
        int target=static_cast<int>(uniform(rng)*available),destination=0;
        while(target>=free_slots[destination])target-=free_slots[destination++];
        target=static_cast<int>(uniform(rng)*residual);
        int state=0;
        while(target>=pending[state])target-=pending[state++];
        --pending[state];++sites[recipients[destination]][state];
        --free_slots[destination];--residual;--available;
    }
    return removal;
}
}
