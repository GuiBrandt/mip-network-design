
#include "network_design/mip/polynomial.hpp"

namespace network_design {
namespace mip {
namespace polynomial {

mip_vars_t::mip_vars_t(GRBModel& model, const instance_t& data)
    : partition(data.graph), partition_used(data.graph.nodeNum()),
      circuit_partition_node(data.graph), circuit_node(data.graph),
      circuit_order(data.graph), circuit_arc(data.graph),
      star_edge(data.graph) {
    static char var_name[32];

    const auto& G = data.graph;
    const int N_PARTITIONS = G.nodeNum();

    // Incidência nas partições
    for (Graph::NodeIt i(G); i != lemon::INVALID; ++i) {
        partition[i].resize(N_PARTITIONS);
        for (int j = 0; j < N_PARTITIONS; j++) {
            std::snprintf(var_name, sizeof(var_name), "partition[%d,%d]",
                          G.id(i), j);
            partition[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, var_name);
            partition[i][j].set(GRB_IntAttr_PoolIgnore, 1);
        }
    }

    // Uso das partições
    for (int j = 0; j < N_PARTITIONS; j++) {
        std::snprintf(var_name, sizeof(var_name), "partition_used[%d]", j);
        partition_used[j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, var_name);
        partition_used[j].set(GRB_IntAttr_BranchPriority, 1000 * j);
        partition_used[j].set(GRB_DoubleAttr_VarHintVal, 0);
        partition_used[j].set(GRB_IntAttr_PoolIgnore, 1);
    }

    // Nós no circuito
    for (Graph::NodeIt i(G); i != lemon::INVALID; ++i) {
        circuit_partition_node[i].resize(N_PARTITIONS);
        for (int j = 0; j < N_PARTITIONS; j++) {
            std::snprintf(var_name, sizeof(var_name),
                          "circuit_partition_node[%d,%d]", G.id(i), j);
            circuit_partition_node[i][j] =
                model.addVar(0.0, 1.0, 0.0, GRB_BINARY, var_name);
        }

        std::snprintf(var_name, sizeof(var_name), "circuit_node[%d]", G.id(i));
        circuit_node[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, var_name);
        circuit_node[i].set(GRB_IntAttr_BranchPriority, 500);
        circuit_node[i].set(GRB_IntAttr_PoolIgnore, 1);
    }

    // Ordem dos nós no circuito
    for (Graph::NodeIt i(G); i != lemon::INVALID; ++i) {
        std::snprintf(var_name, sizeof(var_name), "circuit_order[%d]", G.id(i));
        circuit_order[i] =
            model.addVar(0.0, N_PARTITIONS - 1, 0.0, GRB_INTEGER, var_name);
        circuit_order[i].set(GRB_IntAttr_BranchPriority, -1000);
        circuit_order[i].set(GRB_IntAttr_PoolIgnore, 1);
    }

    // Arcos no circuito
    for (Graph::ArcIt a(G); a != lemon::INVALID; ++a) {
        Graph::Node u = G.source(a), v = G.target(a);
        std::snprintf(var_name, sizeof(var_name), "circuit_arc[%d,%d]", G.id(u),
                      G.id(v));
        uint64_t cost = data.circuit_cost_factor * data.edge_cost[G.edge(u, v)];
        circuit_arc[a] = model.addVar(0.0, 1.0, cost, GRB_BINARY, var_name);
    }

    // Arestas nas estrelas
    for (Graph::EdgeIt e(G); e != lemon::INVALID; ++e) {
        std::snprintf(var_name, sizeof(var_name), "star_edge[%d,%d]",
                      G.id(G.u(e)), G.id(G.v(e)));
        uint64_t cost = data.edge_cost[e];
        star_edge[e] = model.addVar(0.0, 1.0, cost, GRB_BINARY, var_name);
        star_edge[e].set(GRB_IntAttr_BranchPriority, -1000);
    }
}

formulation_t::formulation_t(const instance_t& instance, const GRBEnv& env)
    : instance(instance), model(env), vars(model, instance), G(instance.graph),
      N_PARTITIONS(G.nodeNum()) {
    model.set(GRB_StringAttr_ModelName, "Um Problema de Network Design");
    model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);

    add_used_partition_constraints();
    add_partition_packing_constraints();
    add_partition_symmetry_breaking_constraints();
    add_star_edge_constraints();
    add_circuit_node_constraints();
    add_circuit_arc_constraints();
    add_circuit_order_constraints();
    add_circuit_symmetry_breaking_constraints();

    model.update();
    model.write("model.lp");
}

static char constr_name[64];

void formulation_t::add_used_partition_constraints() {
    for (int j = 0; j < N_PARTITIONS; j++) {
        GRBLinExpr expr;
        for (Graph::NodeIt i(G); i != lemon::INVALID; ++i) {
            std::snprintf(constr_name, sizeof(constr_name),
                          "Partition %d - Used with node %d", j, G.id(i));
            model.addConstr(vars.partition_used[j] >= vars.partition[i][j],
                            constr_name);
        }
    }

    GRBLinExpr expr;
    for (int j = 0; j < N_PARTITIONS; j++) {
        expr += vars.partition_used[j];
    }
    model.addConstr(expr >= 3, "At least 3 partitions");
}

void formulation_t::add_partition_symmetry_breaking_constraints() {
    for (int j = 1; j < N_PARTITIONS; j++) {
        std::snprintf(constr_name, sizeof(constr_name),
                      "Partition %d - Symmetry breaking", j);
        model.addConstr(vars.partition_used[j] <= vars.partition_used[j - 1],
                        constr_name);
    }
}

void formulation_t::add_partition_packing_constraints() {
    // Limite de capacidade.
    for (int j = 0; j < N_PARTITIONS; j++) {
        GRBLinExpr expr;
        for (Graph::NodeIt i(G); i != lemon::INVALID; ++i) {
            expr += instance.node_weight[i] * vars.partition[i][j];
        }
        std::snprintf(constr_name, sizeof(constr_name),
                      "Partition %d - Capacity limit", j);
        model.addConstr(expr <= instance.capacity, constr_name);
    }

    // Um vértice está em exatamente uma partição.
    for (Graph::NodeIt i(G); i != lemon::INVALID; ++i) {
        GRBLinExpr expr;
        for (int j = 0; j < N_PARTITIONS; j++) {
            expr += vars.partition[i][j];
        }
        std::snprintf(constr_name, sizeof(constr_name),
                      "Node %d - Exactly one partition", G.id(i));
        model.addConstr(expr == 1, constr_name);
    }
}

void formulation_t::add_star_edge_constraints() {
    for (int j = 0; j < N_PARTITIONS; j++) {
        for (Graph::EdgeIt e(G); e != lemon::INVALID; ++e) {
            Graph::Node u = G.u(e), v = G.v(e);
            std::snprintf(constr_name, sizeof(constr_name),
                          "Partition %d - Star edge {%d, %d}", j, G.id(u),
                          G.id(v));
            model.addConstr(vars.star_edge[e] >=
                                vars.circuit_partition_node[u][j] -
                                    2 * (1 - vars.partition[v][j]),
                            constr_name);
            model.addConstr(vars.star_edge[e] >=
                                vars.circuit_partition_node[v][j] -
                                    2 * (1 - vars.partition[u][j]),
                            constr_name);
        }
    }
}

void formulation_t::add_circuit_node_constraints() {
    for (int j = 0; j < N_PARTITIONS; j++) {
        // Exatamente um vértice da partição está no circuito se ela for
        // usada, ou nenhum, caso contrário.
        GRBLinExpr nodes_in_partition_expr;
        for (Graph::NodeIt i(G); i != lemon::INVALID; ++i) {
            nodes_in_partition_expr += vars.circuit_partition_node[i][j];
        }
        std::snprintf(constr_name, sizeof(constr_name),
                      "Partition %d - Exactly one circuit node", j);
        model.addConstr(nodes_in_partition_expr == vars.partition_used[j],
                        constr_name);

        // Um vértice só pode estar no circuito em uma partição se de fato
        // estiver naquela partição.
        for (Graph::NodeIt i(G); i != lemon::INVALID; ++i) {
            std::snprintf(constr_name, sizeof(constr_name),
                          "Partition %d - Circuit node %d validity", j,
                          G.id(i));
            model.addConstr(vars.circuit_partition_node[i][j] <=
                                vars.partition[i][j],
                            constr_name);
        }
    }

    for (Graph::NodeIt i(G); i != lemon::INVALID; ++i) {
        for (int j = 0; j < N_PARTITIONS; j++) {
            model.addConstr(vars.circuit_node[i] >=
                            vars.circuit_partition_node[i][j]);
        }
    }
}

void formulation_t::add_circuit_arc_constraints() {
    // Um arco só pode estar no circuito se ambos seus extremos estão no
    // circuito.
    for (Graph::ArcIt a(G); a != lemon::INVALID; ++a) {
        Graph::Node s = G.source(a), t = G.target(a);
        std::snprintf(constr_name, sizeof(constr_name),
                      "Circuit arc %d - Both ends in circuit", G.id(a));
        model.addConstr(2 * vars.circuit_arc[a] <=
                            vars.circuit_node[s] + vars.circuit_node[t],
                        constr_name);
    }

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
                      "Circuit node %d - In-Degree 1", G.id(u));
        model.addConstr(in_degree_expr == vars.circuit_node[u], constr_name);

        std::snprintf(constr_name, sizeof(constr_name),
                      "Circuit node %d - Out-Degree 1", G.id(u));
        model.addConstr(out_degree_expr == vars.circuit_node[u], constr_name);
    }
}

void formulation_t::add_circuit_order_constraints() {
    // Limite da ordem cíclica.
    GRBLinExpr expr;
    for (int j = 0; j < N_PARTITIONS; j++) {
        expr += vars.partition_used[j];
    }
    for (Graph::NodeIt i(G); i != lemon::INVALID; ++i) {
        std::snprintf(constr_name, sizeof(constr_name),
                      "Circuit order of node %d - Number of partitions",
                      G.id(i));
        model.addConstr(vars.circuit_order[i] <= expr - 1);
    }

    // Eliminação de subciclo: se uma aresta está no circuito, então um de
    // seus extremos vem antes do outro na ordem cíclica, exceto para o nó
    // na primeira partição.
    for (Graph::ArcIt a(G); a != lemon::INVALID; ++a) {
        std::snprintf(constr_name, sizeof(constr_name),
                      "Circuit arc %d - Order", G.id(a));
        Graph::Node u = G.source(a), v = G.target(a);
        model.addConstr(vars.circuit_order[v] >=
                            vars.circuit_order[u] + 1 -
                                N_PARTITIONS * (1 - vars.circuit_arc[a]) -
                                N_PARTITIONS *
                                    vars.circuit_partition_node[v][0],
                        constr_name);
    }
}

void formulation_t::add_circuit_symmetry_breaking_constraints() {
    // O vértice na primeira partição tem ID menor que todos os outros
    // vértices no ciclo. Elimina simetrias de rotação.
    for (Graph::EdgeIt e(G); e != lemon::INVALID; ++e) {
        Graph::Node u = G.u(e), v = G.v(e);

        if (G.id(u) > G.id(v)) {
            std::swap(u, v);
        }

        model.addConstr(vars.circuit_partition_node[u][0] >=
                            vars.circuit_partition_node[v][0] -
                                (1 - vars.circuit_node[u]),
                        constr_name);
    }
}

solution_t build_solution(const instance_t& data, const mip_vars_t& vars) {
    solution_t solution(data);

    const auto& G = data.graph;

    // Partições
    int max_partition = 0;
    for (Graph::NodeIt i(G); i != lemon::INVALID; ++i) {
        for (int j = 0; j < G.nodeNum(); j++) {
            if (vars.partition[i][j].get(GRB_DoubleAttr_X) >= 0.5) {
                solution.partition[i] = j;
                if (j > max_partition) {
                    max_partition = j;
                }
            }
        }
    }

    // Nós do circuito
    solution.circuit_nodes.reserve(max_partition + 1);
    for (Graph::NodeIt v(G); v != lemon::INVALID; ++v) {
        for (int j = 0; j <= max_partition; j++) {
            if (vars.circuit_partition_node[v][j].get(GRB_DoubleAttr_X) >=
                0.5) {
                solution.circuit_nodes.push_back(v);
                break;
            }
        }
        if (solution.circuit_nodes.size() > max_partition) {
            break;
        }
    }

    assert(solution.circuit_nodes.size() == max_partition + 1);

    // Ordena os nós de acordo com a ordem do circuito
    std::sort(solution.circuit_nodes.begin(), solution.circuit_nodes.end(),
              [&vars](const Graph::Node& u, const Graph::Node& v) {
                  return vars.circuit_order[u].get(GRB_DoubleAttr_X) <=
                         vars.circuit_order[v].get(GRB_DoubleAttr_X);
              });

    // Índice reverso para a ordem das partições no circuito
    solution.partition_repr.resize(max_partition + 1);
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
