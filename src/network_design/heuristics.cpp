#include "network_design/heuristics.hpp"

namespace network_design {

std::vector<Graph::Node> pack_star(Graph::Node u,
                                   const std::set<Graph::Node>& node_pool,
                                   const instance_t& data) {
    const auto& G = data.graph;

    std::vector<Graph::Node> neighbors;
    for (auto v : node_pool) {
        if (u != v) {
            neighbors.push_back(v);
        }
    }

    if (neighbors.empty()) {
        return {u};
    }

    std::vector<Graph::Node> star = {u};
    std::sort(neighbors.begin(), neighbors.end(),
              [&](const Graph::Node& x, const Graph::Node& y) {
                  return data.edge_cost[G.edge(u, x)] * data.node_weight[x] <
                         data.edge_cost[G.edge(u, y)] * data.node_weight[y];
              });

    uint64_t used = data.node_weight[u];
    for (auto v : neighbors) {
        used += data.node_weight[v];
        if (used > data.capacity) {
            break;
        }
        star.push_back(v);
    }

    return star;
}

solution_t greedy_heuristic(const instance_t& data) {
    const auto& G = data.graph;

    std::set<Graph::Node> node_pool;
    for (Graph::NodeIt v(G); v != lemon::INVALID; ++v) {
        node_pool.insert(v);
    }

    solution_t solution(data);

    int part = -1;
    while (!node_pool.empty()) {
        uint64_t best_cost = std::numeric_limits<uint64_t>::max();
        std::vector<Graph::Node> best_star;
        Graph::Node best_center;

        for (auto u : node_pool) {
            auto star = pack_star(u, node_pool, data);
            uint64_t cost = 0;
            for (auto v : star) {
                if (u != v) {
                    cost += data.edge_cost[G.edge(u, v)];
                }
            }
            if (!solution.circuit_nodes.empty()) {
                const auto prev =
                    solution.circuit_nodes[solution.circuit_nodes.size() - 1];
                cost +=
                    data.circuit_cost_factor * data.edge_cost[G.edge(prev, u)];
            }
            if (cost < best_cost) {
                best_cost = cost;
                best_star = std::move(star);
                best_center = u;
            }
        }

        solution.partition[best_center] = ++part;
        solution.partition_repr.push_back(best_center);
        solution.circuit_nodes.push_back(best_center);
        for (auto v : best_star) {
            node_pool.erase(v);
            solution.partition[v] = part;
        }
    }

    return solution;
}

}; // namespace network_design