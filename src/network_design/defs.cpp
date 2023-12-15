#include "network_design/defs.hpp"

#include <limits>

namespace network_design {

instance_t::instance_t(int nnodes)
    : graph(nnodes), edge_cost(graph), node_weight(graph) {}

solution_t::solution_t(const instance_t& data)
    : instance(data), partition(instance.graph) {}

solution_t::solution_t(const solution_t& other)
    : instance(other.instance), partition(instance.graph),
      circuit_nodes(other.circuit_nodes) {
    // As classes de `Map` do Lemon não implementam construtor de cópia, então
    // é necessário fazer a cópia manualmente.
    for (Graph::NodeIt i(instance.graph); i != lemon::INVALID; ++i) {
        partition[i] = other.partition[i];
    }
}

std::vector<Graph::Arc> solution_t::circuit_arcs() const {
    const size_t N = circuit_nodes.size();
    std::vector<Graph::Arc> result;
    result.reserve(N);
    for (int i = 0; i < N; i++) {
        result.push_back(
            instance.graph.arc(circuit_nodes[i], circuit_nodes[(i + 1) % N]));
    }
    return result;
}

std::vector<Graph::Arc> solution_t::star_arcs() const {
    std::vector<Graph::Arc> result;
    result.reserve(lemon::countNodes(instance.graph) - circuit_nodes.size());
    for (Graph::NodeIt u(instance.graph); u != lemon::INVALID; ++u) {
        const auto& v = circuit_nodes[partition[u]];
        if (u != v) {
            result.push_back(instance.graph.arc(v, u));
        }
    }
    return result;
}

double solution_t::cost() const {
    double value = 0;
    for (const auto& a : circuit_arcs()) {
        value += instance.circuit_cost_factor * instance.edge_cost[a];
    }
    for (const auto& a : star_arcs()) {
        value += instance.edge_cost[a];
    }
    return value;
}

}; // namespace network_design
