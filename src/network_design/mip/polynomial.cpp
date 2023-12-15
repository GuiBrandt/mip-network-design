#include "network_design/mip/polynomial.hpp"

#include <limits>

namespace network_design {
namespace mip {
namespace polynomial {

mip_vars_t::mip_vars_t(GRBModel& model, const instance_t& instance)
    : node_order(instance.graph), mip_vars_base(model, instance) {
    static char var_name[32];

    const auto& G = instance.graph;

    // Ordem dos nós no circuito
    for (Graph::NodeIt i(G); i != lemon::INVALID; ++i) {
        std::snprintf(var_name, sizeof(var_name), "node_order[%d]", G.id(i));
        node_order[i] =
            model.addVar(0.0, G.nodeNum() - 1, 0.0, GRB_CONTINUOUS, var_name);
    }
}

formulation_t::formulation_t(const instance_t& instance, const GRBEnv& env)
    : formulation_base(instance, env) {
    add_circuit_order_constraints();
}

static char constr_name[64];

void formulation_t::add_circuit_order_constraints() {
    const auto& G = instance.graph;

    auto v0 = G.nodeFromId(0);

    // Eliminação de subciclo: se uma aresta está no circuito, então um de
    // seus extremos vem antes do outro na ordem cíclica, exceto para o nó
    // v0 ou o nó que atende v0.
    // Versão fortalecida da restrição Miller-Tucker-Zemlin, por Desrochers &
    // Laporte (1991).
    for (Graph::ArcIt a(G); a != lemon::INVALID; ++a) {
        auto s = G.source(a), t = G.target(a);
        if (t == v0) {
            continue;
        }
        std::snprintf(constr_name, sizeof(constr_name),
                      "node_order_arc[%d][%d]", G.id(s), G.id(t));
        model.addConstr(vars.node_order[s] - vars.node_order[t] +
                            G.nodeNum() * vars.circuit_arc[a] +
                            G.nodeNum() * vars.star_arc[G.arc(t, v0)] +
                            (G.nodeNum() - 2) *
                                vars.circuit_arc[G.oppositeArc(a)] <=
                        G.nodeNum() - 1);
    }
}

} // namespace polynomial
} // namespace mip
} // namespace network_design
