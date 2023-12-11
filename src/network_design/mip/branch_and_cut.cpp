#include "network_design/mip/branch_and_cut.hpp"

#include <limits>

#include <lemon/gomory_hu.h>
#include <lemon/list_graph.h>

namespace network_design {
namespace mip {
namespace branch_and_cut {

static char var_name[32];

mip_vars_t::mip_vars_t(GRBModel& model, const instance_t& data)
    : partition(data.graph), circuit_node(data.graph), circuit_edge(data.graph),
      star_arc(data.graph) {

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

    // Nós no circuito
    for (Graph::NodeIt i(G); i != lemon::INVALID; ++i) {
        std::snprintf(var_name, sizeof(var_name), "circuit_node[%d]", G.id(i));
        circuit_node[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, var_name);
        circuit_node[i].set(GRB_IntAttr_PoolIgnore, 1);
    }

    // Arcos no circuito
    for (Graph::EdgeIt e(G); e != lemon::INVALID; ++e) {
        std::snprintf(var_name, sizeof(var_name), "circuit_edge[%d,%d]",
                      G.id(G.u(e)), G.id(G.v(e)));
        double cost = data.circuit_cost_factor * data.edge_cost[e];

        // Desabilita arestas sem peso definido
        bool usable = cost < 1e99;
        circuit_edge[e] = model.addVar(0.0, usable, cost, GRB_BINARY, var_name);
    }

    // Arcos nas estrelas
    for (Graph::ArcIt a(G); a != lemon::INVALID; ++a) {
        std::snprintf(var_name, sizeof(var_name), "star_arc[%d,%d]",
                      G.id(G.source(a)), G.id(G.target(a)));
        double cost = data.edge_cost[a];

        // Desabilita arestas sem peso definido
        bool usable = cost < 1e99;
        star_arc[a] = model.addVar(0.0, usable, cost, GRB_BINARY, var_name);
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
    add_star_arc_constraints();
    add_circuit_edge_constraints();

    model.update();
    model.write("model.lp");
}

static char constr_name[64];

void formulation_t::add_used_partition_constraints() {
    for (Graph::NodeIt u(G); u != lemon::INVALID; ++u) {
        GRBLinExpr expr;
        for (Graph::NodeIt v(G); v != lemon::INVALID; ++v) {
            std::snprintf(constr_name, sizeof(constr_name),
                          "Partition %d - Used with node %d", G.id(u), G.id(v));
            model.addConstr(vars.circuit_node[u] >= vars.partition[v][G.id(u)],
                            constr_name);
        }
        model.addConstr(vars.partition[u][G.id(u)] >= vars.circuit_node[u]);
    }

    GRBLinExpr expr;
    for (Graph::NodeIt u(G); u != lemon::INVALID; ++u) {
        expr += vars.circuit_node[u];
    }
    model.addConstr(expr >= 3, "At least 3 partitions");
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

void formulation_t::add_star_arc_constraints() {
    for (Graph::ArcIt a(G); a != lemon::INVALID; ++a) {
        Graph::Node s = G.source(a), t = G.target(a);
        std::snprintf(constr_name, sizeof(constr_name), "Star arc (%d, %d)",
                      G.id(s), G.id(t));
        model.addConstr(vars.star_arc[a] >=
                            vars.circuit_node[t] -
                                2 * (1 - vars.partition[s][G.id(t)]),
                        constr_name);
    }
    for (Graph::NodeIt v(G); v != lemon::INVALID; ++v) {
        std::snprintf(constr_name, sizeof(constr_name),
                      "Star node %d out-degree 1", G.id(v));
        GRBLinExpr expr;
        for (Graph::OutArcIt a(G, v); a != lemon::INVALID; ++a) {
            expr += vars.star_arc[a];
        }
        model.addConstr(expr == 1 - vars.circuit_node[v], constr_name);
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
        for (Graph::IncEdgeIt e(G, u); e != lemon::INVALID; ++e) {
            degree_expr += vars.circuit_edge[e];
        }

        std::snprintf(constr_name, sizeof(constr_name),
                      "Circuit node %d - Degree 2", G.id(u));
        model.addConstr(degree_expr == 2 * vars.circuit_node[u], constr_name);
    }
}

void formulation_t::callback() {
    if (where != GRB_CB_MIPSOL) {
        return;
    }

    Graph::EdgeMap<double> capacity(G);
    for (Graph::EdgeIt e(G); e != lemon::INVALID; ++e) {
        Graph::Node u = G.u(e), v = G.v(e);
        capacity[e] = getSolution(vars.circuit_edge[G.edge(u, v)]) +
                      getSolution(vars.star_arc[G.arc(u, v)]) +
                      getSolution(vars.star_arc[G.arc(v, u)]);
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
            Graph::Node u = G.u(e), v = G.v(e);
            if (cutmap[u] != cutmap[v]) {
                expr += vars.circuit_edge[G.edge(u, v)] +
                        vars.star_arc[G.arc(u, v)] + vars.star_arc[G.arc(v, u)];
            }
        }
        addLazy(expr >= 1);
    }
}

solution_t build_solution(const instance_t& data, const mip_vars_t& vars) {
    solution_t solution(data);

    const auto& G = data.graph;

    std::cout << "digraph D {" << std::endl;
    for (Graph::ArcIt a(G); a != lemon::INVALID; ++a) {
        Graph::Node s = G.source(a), t = G.target(a);
        if (vars.star_arc[a].get(GRB_DoubleAttr_X) >= 0.5) {
            std::cout << "\t" << G.id(s) << " -> " << G.id(t)
                      << " [color=blue];" << std::endl;
        }
    }
    for (Graph::EdgeIt e(G); e != lemon::INVALID; ++e) {
        Graph::Node s = G.u(e), t = G.v(e);
        if (vars.circuit_edge[e].get(GRB_DoubleAttr_X) >= 0.5) {
            std::cout << "\t" << G.id(s) << " -> " << G.id(t)
                      << " [color=red, dir=none];" << std::endl;
        }
    }
    std::cout << "}" << std::endl;

    // Partições
    std::vector<int> partition_map(G.nodeNum());
    for (int j = 0; j < G.nodeNum(); j++) {
        partition_map[j] = j;
    }

    int max_part = -1;
    for (Graph::NodeIt v(G); v != lemon::INVALID; ++v) {
        if (vars.circuit_node[v].get(GRB_DoubleAttr_X) >= 0.5) {
            partition_map[G.id(v)] = ++max_part;
        }
    }

    for (Graph::NodeIt i(G); i != lemon::INVALID; ++i) {
        for (int j = 0; j < G.nodeNum(); j++) {
            if (vars.partition[i][j].get(GRB_DoubleAttr_X) >= 0.5) {
                solution.partition[i] = partition_map[j];
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
                vars.circuit_edge[a].get(GRB_DoubleAttr_X) >= 0.5) {
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

} // namespace branch_and_cut
} // namespace mip
} // namespace network_design
