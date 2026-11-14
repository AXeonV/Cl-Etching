#pragma once
#include <map>
#include <string>
namespace etch::app {
using Arguments=std::map<std::string,std::string>;
Arguments parse(int argc,char** argv);
}
