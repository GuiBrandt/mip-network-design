#ifndef _NETWORK_DESIGN_UTILS_HPP
#define _NETWORK_DESIGN_UTILS_HPP

#include "defs.hpp"

#include <fstream>
#include <random>

namespace network_design {

void view(const solution_t&);

std::unique_ptr<instance_t> random_instance(const Graph&, int seed);

std::unique_ptr<instance_t> read_instance(std::string filename);

}; // namespace network_design

#endif // _NETWORK_DESIGN_UTILS_HPP
