#include "network_design/defs.hpp"

namespace network_design {

instance_t::instance_t(const Graph& graph)
    : graph(graph), edge_cost(graph), node_weight(graph) {}

instance_t::instance_t(const instance_t& other)
    : graph(other.graph), edge_cost(graph), node_weight(graph),
      capacity(other.capacity), circuit_cost_factor(other.circuit_cost_factor) {
    // As classes de `Map` do Lemon não implementam construtor de cópia, então
    // é necessário fazer a cópia manualmente.
    for (Graph::EdgeIt e(graph); e != lemon::INVALID; ++e) {
        edge_cost[e] = other.edge_cost[e];
    }
    for (Graph::NodeIt v(graph); v != lemon::INVALID; ++v) {
        node_weight[v] = other.node_weight[v];
    }
}

solution_t::solution_t(const instance_t& data)
    : instance(data), partition(instance.graph) {}

solution_t::solution_t(const solution_t& other)
    : instance(other.instance), partition_repr(other.partition_repr),
      partition(instance.graph), circuit_nodes(other.circuit_nodes) {
    // As classes de `Map` do Lemon não implementam construtor de cópia, então
    // é necessário fazer a cópia manualmente.
    for (Graph::NodeIt i(instance.graph); i != lemon::INVALID; ++i) {
        partition[i] = other.partition[i];
    }
}

solution_t::solution_t(const solution_t&& other)
    : instance(other.instance), partition_repr(std::move(other.partition_repr)),
      partition(instance.graph), circuit_nodes(std::move(other.circuit_nodes)) {
    // As classes de `Map` do Lemon não implementam construtor de cópia, então
    // é necessário fazer a cópia manualmente.
    for (Graph::NodeIt i(instance.graph); i != lemon::INVALID; ++i) {
        partition[i] = other.partition[i];
    }
}

std::vector<Graph::Edge> solution_t::circuit_edges() const {
    const size_t N = circuit_nodes.size();
    std::vector<Graph::Edge> result;
    result.reserve(N);
    for (int i = 0; i < N; i++) {
        result.push_back(
            instance.graph.edge(circuit_nodes[i], circuit_nodes[(i + 1) % N]));
    }
    return result;
}

std::vector<Graph::Edge> solution_t::star_edges() const {
    std::vector<Graph::Edge> result;
    result.reserve(instance.graph.nodeNum() - circuit_nodes.size());
    for (Graph::NodeIt u(instance.graph); u != lemon::INVALID; ++u) {
        const auto& v = partition_repr[partition[u]];
        if (u != v) {
            result.push_back(instance.graph.edge(u, v));
        }
    }
    return result;
}

uint64_t solution_t::cost() const {
    uint64_t value = 0;
    for (const auto& e : circuit_edges()) {
        value += instance.circuit_cost_factor * instance.edge_cost[e];
    }
    for (const auto& e : star_edges()) {
        value += instance.edge_cost[e];
    }
    return value;
}

}; // namespace network_design
