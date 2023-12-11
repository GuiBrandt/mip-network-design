#include "network_design/mip/branch_and_cut.hpp"

#include <limits>

#include <lemon/gomory_hu.h>

namespace network_design {
namespace mip {
namespace branch_and_cut {

mip_vars_t::mip_vars_t(GRBModel& model, const instance_t& data)
    : partition(data.graph), partition_used(data.graph.nodeNum()),
      circuit_partition_node(data.graph), circuit_node(data.graph),
      circuit_edge(data.graph), star_edge(data.graph) {
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
            circuit_partition_node[i][j].set(GRB_IntAttr_PoolIgnore, 1);
        }

        std::snprintf(var_name, sizeof(var_name), "circuit_node[%d]", G.id(i));
        circuit_node[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, var_name);
        circuit_node[i].set(GRB_IntAttr_PoolIgnore, 1);
    }

    // Arestas no circuito
    for (Graph::EdgeIt e(G); e != lemon::INVALID; ++e) {
        Graph::Node u = G.u(e), v = G.v(e);
        std::snprintf(var_name, sizeof(var_name), "circuit_edge[%d,%d]",
                      G.id(u), G.id(v));
        double cost = data.circuit_cost_factor * data.edge_cost[G.edge(u, v)];

        // Desabilita arestas sem peso definido
        bool usable = cost < 1e99;
        circuit_edge[e] = model.addVar(0.0, usable, cost, GRB_BINARY, var_name);
    }

    // Arestas nas estrelas
    for (Graph::EdgeIt e(G); e != lemon::INVALID; ++e) {
        std::snprintf(var_name, sizeof(var_name), "star_edge[%d,%d]",
                      G.id(G.u(e)), G.id(G.v(e)));
        double cost = data.edge_cost[e];

        // Desabilita arestas sem peso definido
        bool usable = cost < 1e99;
        star_edge[e] = model.addVar(0.0, usable, cost, GRB_BINARY, var_name);
        star_edge[e].set(GRB_IntAttr_BranchPriority, -1000);
    }
}

formulation_t::formulation_t(const instance_t& instance, const GRBEnv& env)
    : instance(instance), model(env), vars(model, instance), G(instance.graph),
      N_PARTITIONS(G.nodeNum()) {
    model.set(GRB_StringAttr_ModelName, "Um Problema de Network Design");
    model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
    model.set(GRB_IntParam_LazyConstraints, 1);
    model.setCallback(this);

    add_used_partition_constraints();
    add_partition_packing_constraints();
    add_partition_symmetry_breaking_constraints();
    add_star_edge_constraints();
    add_circuit_node_constraints();
    add_circuit_edge_constraints();
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

void formulation_t::add_circuit_edge_constraints() {
    // Uma aresta só pode estar no circuito se ambos seus extremos estão no
    // circuito.
    for (Graph::EdgeIt e(G); e != lemon::INVALID; ++e) {
        Graph::Node s = G.u(e), t = G.v(e);
        std::snprintf(constr_name, sizeof(constr_name),
                      "Circuit edge %d - Both ends in circuit", G.id(e));
        model.addConstr(2 * vars.circuit_edge[e] <=
                            vars.circuit_node[s] + vars.circuit_node[t],
                        constr_name);
    }

    // Todo vértice no circuito deve ter grau 2.
    for (Graph::NodeIt u(G); u != lemon::INVALID; ++u) {
        GRBLinExpr degree_expr;
        for (Graph::NodeIt v(G); v != lemon::INVALID; ++v) {
            if (u != v) {
                degree_expr += vars.circuit_edge[G.edge(u, v)];
            }
        }

        std::snprintf(constr_name, sizeof(constr_name),
                      "Circuit node %d - Degree 2", G.id(u));
        model.addConstr(degree_expr == 2 * vars.circuit_node[u], constr_name);
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

void formulation_t::callback() {
    Graph::EdgeMap<double> capacity(G);

    switch (where) {
    case GRB_CB_MIPSOL: {
        for (Graph::EdgeIt e(G); e != lemon::INVALID; ++e) {
            capacity[e] = getSolution(vars.circuit_edge[e]) +
                          getSolution(vars.star_edge[e]);
        }
        break;
    }
    case GRB_CB_MIPNODE: {
        if (getIntInfo(GRB_CB_MIPNODE_STATUS) != GRB_OPTIMAL) {
            return;
        }
        for (Graph::EdgeIt e(G); e != lemon::INVALID; ++e) {
            capacity[e] = getNodeRel(vars.circuit_edge[e]) +
                          getNodeRel(vars.star_edge[e]);
        }
        break;
    }
    default:
        return;
    }

    lemon::GomoryHu<Graph, Graph::EdgeMap<double>> gomory_hu(G, capacity);
    gomory_hu.run();

    Graph::NodeMap<bool> cutmap(G);
    for (Graph::NodeIt u(G); u != lemon::INVALID; ++u) {
        if (gomory_hu.predNode(u) == lemon::INVALID ||
            gomory_hu.predValue(u) > 1.0 - 1e-5) {
            continue;
        }

        gomory_hu.minCutMap(u, gomory_hu.predNode(u), cutmap);

        GRBLinExpr expr;
        for (Graph::EdgeIt e(G); e != lemon::INVALID; ++e) {
            if (cutmap[G.u(e)] != cutmap[G.v(e)]) {
                expr += vars.circuit_edge[e] + vars.star_edge[e];
            }
        }

        addLazy(expr >= 1);
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
    // TODO:
    // std::sort(solution.circuit_nodes.begin(), solution.circuit_nodes.end(),
    //           [&vars](const Graph::Node& u, const Graph::Node& v) {
    //               return vars.circuit_order[u].get(GRB_DoubleAttr_X) <=
    //                      vars.circuit_order[v].get(GRB_DoubleAttr_X);
    //           });

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

} // namespace branch_and_cut
} // namespace mip
} // namespace network_design
