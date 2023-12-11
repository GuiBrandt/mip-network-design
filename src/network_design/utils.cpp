#include "network_design/utils.hpp"

namespace network_design {

void view(const solution_t& solution) {
    const auto& G = solution.instance.graph;

    std::ofstream out("solution.dot");
    out << "graph {" << std::endl;
    for (auto e : solution.circuit_edges()) {
        out << "\t" << G.id(G.u(e)) << " -- " << G.id(G.v(e))
            << " [ color=red ];" << std::endl;
    }
    for (auto e : solution.star_edges()) {
        out << "\t" << G.id(G.u(e)) << " -- " << G.id(G.v(e))
            << " [ color=blue ];" << std::endl;
    }
    out << "}" << std::endl;
    out.close();

    system("circo -q -Tpdf solution.dot -o solution.pdf");
    system("okular solution.pdf");
}

instance_t random_instance(const Graph& G, int seed) {
    std::minstd_rand rng(seed);

    instance_t instance(G);

    Graph::NodeMap<std::pair<int, int>> position(G);

    std::uniform_int_distribution<uint64_t> capacity_dist(20, 50);
    instance.capacity = capacity_dist(rng);

    std::uniform_int_distribution<uint64_t> factor_dist(2, 5);
    instance.circuit_cost_factor = factor_dist(rng);

    std::uniform_int_distribution<int> position_dist(0, 10);
    std::uniform_int_distribution<uint64_t> weight_dist(1, 10);
    for (Graph::NodeIt v(G); v != lemon::INVALID; ++v) {
        instance.node_weight[v] = weight_dist(rng);
        position[v] = {position_dist(rng), position_dist(rng)};
    }

    for (Graph::EdgeIt e(G); e != lemon::INVALID; ++e) {
        Graph::Node u = G.u(e), v = G.v(e);
        instance.edge_cost[e] = abs(position[u].first - position[v].first) +
                                abs(position[u].second - position[v].second);
    }

    return instance;
}
}; // namespace network_design
