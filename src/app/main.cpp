#include "options.hpp"
#include "planar/backend.hpp"
#include "spatial/backend.hpp"
#include <chrono>
#include <iostream>
#include <stdexcept>
namespace etch::app {
Distribution energy(const Arguments& args) {
    std::map<std::string,std::function<Distribution()>> factories={
        {"constant",[&]{return constant(std::stod(args.at("--energy-ev")));}},
        {"gaussian",[&]{return gaussian(std::stod(args.at("--energy-ev")),std::stod(args.at("--energy-sigma")),
                                        0,std::numeric_limits<double>::infinity());}}
    };
    return factories.at(args.at("--energy"))();
}
std::unique_ptr<Reaction> reaction(Arguments& args) {
    if(args.at("--model")=="chlorination") { args["--angular"]="not_used";return std::make_unique<ChlorinationReaction>(); }
    if(args.at("--model")!="yield")throw std::invalid_argument("Unknown reaction model");
    std::map<std::string,AngularYield> laws={{"cosine",cosine_yield},{"enhanced",enhanced_yield}};
    return std::make_unique<YieldReaction>(laws.at(args.at("--angular")));
}
DirectionDistribution direction(const Arguments& args) {
    std::map<std::string,std::function<DirectionDistribution()>> factories={
        {"gaussian",[&]{return gaussian_direction(std::stod(args.at("--fwhm"))*std::acos(-1.0)/180/2.354820045);}},
        {"cosine",[]{return cosine_direction();}},{"normal",[]{return normal_direction();}}
    };
    return factories.at(args.at("--angle"))();
}

template<class Backend> void run(typename Backend::Geometry geometry,typename Backend::Source source,Arguments& args) {
    geometry.cell_nm=std::stod(args.at("--cell-nm"));
    RunSettings<typename Backend::Geometry> settings{geometry};
    settings.seed=std::stoull(args.at("--seed"));settings.reflections=std::stoi(args.at("--reflect"));
    settings.surface_subdivisions=std::stoi(args.at("--surface-subdivisions"));
    settings.max_bounces=std::stoi(args.at("--max-bounces"));
    Engine<Backend> simulation(settings,source,reaction(args));
    auto start=std::chrono::steady_clock::now();
    long long particles=std::stoll(args.at("--particles")),batch=std::min(250000LL,particles);
    for(long long done=0;done<particles;done+=batch) {
        simulation.advance(std::min(batch,particles-done));
        std::cout<<simulation.stats.particles<<" particles, "<<simulation.reaction->removed
                 <<" Si removed, "<<simulation.stats.deleted_cells<<" cells"<<std::endl;
    }
    double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    simulation.save(args.at("--output"),elapsed);
    std::ofstream config(args.at("--output")+"/arguments.txt");
    for(const auto& [key,value]:args)config<<key<<' '<<value<<'\n';
    std::cout<<"Saved "<<args.at("--output")<<" ("<<elapsed<<" s)"<<std::endl;
}
}
int main(int argc,char** argv) {
    using namespace etch;using namespace etch::app;
    auto args=parse(argc,argv);int dimension=std::stoi(args.at("--dimension"));
    int nx=std::stoi(args.at("--nx")),nz=std::stoi(args.at("--nz")),z=std::stoi(args.at("--surface-z"));
    int mask=std::stoi(args.at("--mask-thickness")),opening=std::stoi(args.at("--opening"));
    if(dimension==2)run<planar::Backend>({nx,nz,z,mask,opening},planar::Source{direction(args),energy(args)},args);
    else if(dimension==3) {
        int ny=std::stoi(args.at("--ny"));
        int opening_y=args.at("--opening-y").empty()?ny:std::stoi(args.at("--opening-y"));
        run<spatial::Backend>({nx,ny,nz,z,mask,opening,opening_y},spatial::Source{direction(args),energy(args)},args);
    } else throw std::invalid_argument("--dimension must be 2 or 3");
}
