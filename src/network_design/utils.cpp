#include "network_design/utils.hpp"

#include <mylib/mygraphlib.h>

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

std::unique_ptr<instance_t> random_instance(int nnodes, int seed) {
    std::minstd_rand rng(seed);

    std::unique_ptr<instance_t> instance(new instance_t(nnodes));

    Graph::NodeMap<std::pair<int, int>> position(instance->graph);

    std::uniform_real_distribution<double> capacity_dist(20, 50);
    instance->capacity = capacity_dist(rng);

    std::uniform_real_distribution<double> factor_dist(2, 5);
    instance->circuit_cost_factor = factor_dist(rng);

    std::uniform_real_distribution<double> position_dist(0, 10);
    std::uniform_real_distribution<double> weight_dist(1, 10);
    for (Graph::NodeIt v(instance->graph); v != lemon::INVALID; ++v) {
        instance->node_weight[v] = weight_dist(rng);
        position[v] = {position_dist(rng), position_dist(rng)};
    }

    for (Graph::EdgeIt e(instance->graph); e != lemon::INVALID; ++e) {
        Graph::Node u = instance->graph.u(e), v = instance->graph.v(e);
        instance->edge_cost[e] = abs(position[u].first - position[v].first) +
                                 abs(position[u].second - position[v].second);
    }

    return instance;
}

std::unique_ptr<instance_t> read_instance(std::string filename) {
    using LGraph = lemon::ListGraph;

    LGraph H;
    GraphTable GT(filename, H);

    int nnodes;
    GT.Get("nnodes", nnodes);

    std::unique_ptr<instance_t> instance(new instance_t(nnodes));

    GT.Get("capacity", instance->capacity);
    GT.Get("gamma", instance->circuit_cost_factor);

    LGraph::NodeMap<double> weight(H);
    GT.Get("weight", weight);
    for (LGraph::NodeIt v(H); v != lemon::INVALID; ++v) {
        instance->node_weight[instance->graph.nodeFromId(H.id(v))] = weight[v];
    }

    LGraph::EdgeMap<double> cost(H);
    GT.Get("distance", cost);
    for (LGraph::EdgeIt e(H); e != lemon::INVALID; ++e) {
        auto u = instance->graph.nodeFromId(H.id(H.u(e))),
             v = instance->graph.nodeFromId(H.id(H.v(e)));
        instance->edge_cost[instance->graph.edge(u, v)] = cost[e];
    }

    return instance;
}

}; // namespace network_design
