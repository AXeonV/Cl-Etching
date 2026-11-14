#include "options.hpp"
namespace etch::app {
Arguments parse(int argc,char** argv) {
    Arguments args={
        {"--dimension","2"},{"--model","yield"},{"--angular","cosine"},{"--angle","gaussian"},
        {"--energy","constant"},{"--seed","42"},{"--output","outputs/run"},
        {"--reflect","1"},{"--max-bounces","0"},{"--nx",""},{"--ny",""},{"--nz",""},
        {"--surface-z",""},{"--opening",""},{"--opening-y",""},{"--mask-thickness",""},
        {"--particles",""},{"--fwhm",""},{"--energy-ev",""},{"--energy-sigma",""}
    };
    for(int i=1;i<argc;i+=2)args.at(argv[i])=argv[i+1];
    return args;
}
}
