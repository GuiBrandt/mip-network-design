#include "network_design/mip/base.hpp"

namespace network_design {
namespace mip {

mip_vars_base::mip_vars_base(GRBModel& model, const instance_t& data)
    : circuit_node(data.graph), circuit_arc(data.graph), star_arc(data.graph) {
    static char var_name[32];

    const auto& G = data.graph;
    const int N_PARTITIONS = G.nodeNum();

    // Nós no circuito
    for (Graph::NodeIt i(G); i != lemon::INVALID; ++i) {
        std::snprintf(var_name, sizeof(var_name), "circuit_node[%d]", G.id(i));
        circuit_node[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, var_name);
    }

    // Arcos no circuito
    for (Graph::ArcIt a(G); a != lemon::INVALID; ++a) {
        Graph::Node u = G.source(a), v = G.target(a);
        std::snprintf(var_name, sizeof(var_name), "circuit_arc[%d,%d]", G.id(u),
                      G.id(v));
        double cost = data.circuit_cost_factor * data.edge_cost[G.edge(u, v)];

        // Desabilita arestas sem peso definido
        bool usable = cost > 0;
        circuit_arc[a] = model.addVar(0.0, usable, cost, GRB_BINARY, var_name);
    }

    // Arestas nas estrelas
    for (Graph::ArcIt a(G); a != lemon::INVALID; ++a) {
        std::snprintf(var_name, sizeof(var_name), "star_arc[%d,%d]",
                      G.id(G.source(a)), G.id(G.target(a)));
        double cost = data.edge_cost[a];

        // Desabilita arestas sem peso definido
        bool usable = cost < 1e99;
        star_arc[a] = model.addVar(0.0, usable, cost, GRB_BINARY, var_name);
    }
}

} // namespace mip
} // namespace network_design
