#ifndef _NETWORK_DESIGN_HEURISTICS_HPP
#define _NETWORK_DESIGN_HEURISTICS_HPP

#include <limits>
#include <set>
#include <vector>

#include "defs.hpp"

namespace network_design {

/**
 * Heurística gulosa para o problema.
 *
 * Prioriza maximizar o empacotamento das estrelas com o menor custo possível,
 * e tenta minimizar o custo das arestas no ciclo.
 */
solution_t greedy_heuristic(const instance_t& data);

}; // namespace network_design

#endif // _NETWORK_DESIGN_HEURISTICS_HPP
