#include "network_design/mip/polynomial.hpp"

#include <limits>

namespace network_design {
namespace mip {
namespace polynomial {

mip_vars_t::mip_vars_t(GRBModel& model, const instance_t& data)
    : circuit_node(data.graph), circuit_source(data.graph),
      circuit_order(data.graph), circuit_arc(data.graph), star_arc(data.graph) {
    static char var_name[32];

    const auto& G = data.graph;
    const int N_PARTITIONS = G.nodeNum();

    // Nós no circuito
    for (Graph::NodeIt i(G); i != lemon::INVALID; ++i) {
        std::snprintf(var_name, sizeof(var_name), "circuit_node[%d]", G.id(i));
        circuit_node[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, var_name);
        circuit_node[i].set(GRB_IntAttr_PoolIgnore, 1);
    }

    for (Graph::NodeIt i(G); i != lemon::INVALID; ++i) {
        std::snprintf(var_name, sizeof(var_name), "circuit_source[%d]",
                      G.id(i));
        circuit_source[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, var_name);
        circuit_source[i].set(GRB_IntAttr_PoolIgnore, 1);
    }

    // Ordem dos nós no circuito
    for (Graph::NodeIt i(G); i != lemon::INVALID; ++i) {
        std::snprintf(var_name, sizeof(var_name), "circuit_order[%d]", G.id(i));
        circuit_order[i] =
            model.addVar(0.0, N_PARTITIONS - 1, 0.0, GRB_INTEGER, var_name);
        // circuit_order[i].set(GRB_IntAttr_BranchPriority, -1000);
        circuit_order[i].set(GRB_IntAttr_PoolIgnore, 1);
    }

    // Arcos no circuito
    for (Graph::ArcIt a(G); a != lemon::INVALID; ++a) {
        Graph::Node u = G.source(a), v = G.target(a);
        std::snprintf(var_name, sizeof(var_name), "circuit_arc[%d,%d]", G.id(u),
                      G.id(v));
        double cost = data.circuit_cost_factor * data.edge_cost[G.edge(u, v)];

        // Desabilita arestas sem peso definido
        bool usable = cost < 1e99;
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
        // star_arc[a].set(GRB_IntAttr_BranchPriority, -1000);
    }
}

formulation_t::formulation_t(const instance_t& instance, const GRBEnv& env)
    : instance(instance), model(env), vars(model, instance), G(instance.graph),
      N_PARTITIONS(G.nodeNum()) {
    model.set(GRB_StringAttr_ModelName, "Um Problema de Network Design");
    model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);

    add_partition_constraints();
    add_star_arc_constraints();
    add_circuit_arc_constraints();
    add_circuit_order_constraints();

    model.update();
    model.write("model.lp");
}

static char constr_name[64];

void formulation_t::add_partition_constraints() {
    // Limite de capacidade.
    for (Graph::NodeIt u(G); u != lemon::INVALID; ++u) {
        GRBLinExpr expr = instance.node_weight[u];
        for (Graph::OutArcIt a(G, u); a != lemon::INVALID; ++a) {
            expr += instance.node_weight[G.target(a)] * vars.star_arc[a];
        }
        std::snprintf(constr_name, sizeof(constr_name),
                      "partition_capacity[%d]", G.id(u));
        model.addConstr(expr <= instance.capacity, constr_name);
    }

    // Pelo menos 3 partes são usadas.
    GRBLinExpr expr;
    for (Graph::NodeIt u(G); u != lemon::INVALID; ++u) {
        expr += vars.circuit_node[u];
    }
    model.addConstr(expr >= 3, "At least 3 parts");
}

void formulation_t::add_star_arc_constraints() {
    // Uma aresta de estrela sai de um nó do circuito para um nó fora do
    // circuito.
    for (Graph::ArcIt a(G); a != lemon::INVALID; ++a) {
        Graph::Node s = G.source(a), t = G.target(a);
        model.addConstr(vars.star_arc[a] <= vars.circuit_node[s]);
        model.addConstr(vars.star_arc[a] <= 1 - vars.circuit_node[t]);
    }

    // Exatamente um arco de estrela entra num nó que não está no circuito.
    for (Graph::NodeIt v(G); v != lemon::INVALID; ++v) {
        GRBLinExpr in_degree_expr;
        for (Graph::InArcIt a(G, v); a != lemon::INVALID; ++a) {
            in_degree_expr += vars.star_arc[a];
        }
        std::snprintf(constr_name, sizeof(constr_name),
                      "star_node_in_degree[%d]", G.id(v));
        model.addConstr(in_degree_expr == 1 - vars.circuit_node[v],
                        constr_name);
    }
}

void formulation_t::add_circuit_arc_constraints() {
    // No máximo um dos arcos correspondentes a uma aresta pode ser
    // escolhido para o ciclo.
    for (Graph::NodeIt u(G); u != lemon::INVALID; ++u) {
        for (Graph::NodeIt v(G); v != lemon::INVALID; ++v) {
            if (u != v) {
                model.addConstr(vars.circuit_arc[G.arc(u, v)] +
                                    vars.circuit_arc[G.arc(v, u)] <=
                                1);
            }
        }
    }

    // Todo vértice no circuito deve ter grau de entrada e saída igual a 1.
    for (Graph::NodeIt u(G); u != lemon::INVALID; ++u) {
        GRBLinExpr in_degree_expr;
        GRBLinExpr out_degree_expr;
        for (Graph::NodeIt v(G); v != lemon::INVALID; ++v) {
            if (u != v) {
                in_degree_expr += vars.circuit_arc[G.arc(v, u)];
                out_degree_expr += vars.circuit_arc[G.arc(u, v)];
            }
        }

        std::snprintf(constr_name, sizeof(constr_name),
                      "circuit_node_in_degree[%d]", G.id(u));
        model.addConstr(in_degree_expr == vars.circuit_node[u], constr_name);

        std::snprintf(constr_name, sizeof(constr_name),
                      "circuit_node_out_degree[%d]", G.id(u));
        model.addConstr(out_degree_expr == vars.circuit_node[u], constr_name);
    }
}

void formulation_t::add_circuit_order_constraints() {
    // Fonte do circuito
    GRBLinExpr source_expr;
    for (Graph::NodeIt u(G); u != lemon::INVALID; ++u) {
        source_expr += vars.circuit_source[u];
    }
    std::snprintf(constr_name, sizeof(constr_name), "circuit_source");
    model.addConstr(source_expr == 1);

    for (Graph::NodeIt u(G); u != lemon::INVALID; ++u) {
        model.addConstr(vars.circuit_source[u] <= vars.circuit_node[u]);
    }

    for (Graph::EdgeIt e(G); e != lemon::INVALID; ++e) {
        auto u = G.u(e), v = G.v(e);
        if (G.id(u) < G.id(v)) {
            std::swap(u, v);
        }
        model.addConstr(vars.circuit_source[u] >=
                        vars.circuit_source[v] - (1 - vars.circuit_node[u]));
    }

    // Limite da ordem cíclica.
    GRBLinExpr expr;
    for (Graph::NodeIt u(G); u != lemon::INVALID; ++u) {
        expr += vars.circuit_node[u];
    }
    for (Graph::NodeIt i(G); i != lemon::INVALID; ++i) {
        std::snprintf(constr_name, sizeof(constr_name),
                      "circuit_order_limit[%d]", G.id(i));
        model.addConstr(vars.circuit_order[i] <= expr - 1);
    }

    // Eliminação de subciclo: se uma aresta está no circuito, então um de
    // seus extremos vem antes do outro na ordem cíclica, exceto para o nó
    // na primeira partição.
    for (Graph::ArcIt a(G); a != lemon::INVALID; ++a) {
        Graph::Node s = G.source(a), t = G.target(a);
        std::snprintf(constr_name, sizeof(constr_name),
                      "circuit_order_arc[%d][%d]", G.id(s), G.id(t));
        model.addConstr(vars.circuit_order[t] >=
                            vars.circuit_order[s] + 1 -
                                N_PARTITIONS * (1 - vars.circuit_arc[a]) -
                                N_PARTITIONS * vars.circuit_source[t],
                        constr_name);
    }
}

solution_t build_solution(const instance_t& data, const mip_vars_t& vars) {
    solution_t solution(data);

    const auto& G = data.graph;

    // Partições
    std::vector<int> partition_map(G.nodeNum());
    for (int j = 0; j < G.nodeNum(); j++) {
        partition_map[j] = j;
    }

    int max_part = -1;
    for (Graph::NodeIt v(G); v != lemon::INVALID; ++v) {
        if (vars.circuit_node[v].get(GRB_DoubleAttr_X) >= 0.5) {
            partition_map[G.id(v)] = solution.partition[v] = ++max_part;
        }
    }

    for (Graph::NodeIt u(G); u != lemon::INVALID; ++u) {
        for (Graph::InArcIt a(G, u); a != lemon::INVALID; ++a) {
            if (vars.star_arc[a].get(GRB_DoubleAttr_X) >= 0.5) {
                solution.partition[u] = partition_map[G.id(G.source(a))];
                break;
            }
        }
    }

    // Nós do circuito de acordo com a ordem do circuito
    Graph::Node current, prev;
    for (Graph::NodeIt v(G); v != lemon::INVALID; ++v) {
        if (vars.circuit_node[v].get(GRB_DoubleAttr_X) >= 0.5) {
            current = v;
            break;
        }
    }

    solution.circuit_nodes.reserve(max_part + 1);
    do {
        solution.circuit_nodes.push_back(current);
        for (Graph::OutArcIt a(G, current); a != lemon::INVALID; ++a) {
            auto target = G.target(a);
            if (target != prev &&
                vars.circuit_arc[a].get(GRB_DoubleAttr_X) >= 0.5) {
                prev = current;
                current = target;
                break;
            }
        }
    } while (current != solution.circuit_nodes[0]);

    assert(solution.circuit_nodes.size() == max_part + 1);

    // Índice reverso para a ordem das partições no circuito
    solution.partition_repr.resize(max_part + 1);
    for (auto v : solution.circuit_nodes) {
        solution.partition_repr[solution.partition[v]] = v;
    }

    return solution;
}

solution_t formulation_t::solve() {
    model.optimize();
    model.write("model.sol");
    return build_solution(instance, vars);
}

} // namespace polynomial
} // namespace mip
} // namespace network_design
