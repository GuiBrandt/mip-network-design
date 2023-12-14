#include "network_design/mip/branch_and_cut.hpp"

#include <limits>

#include <lemon/gomory_hu.h>
#include <lemon/list_graph.h>
#include <lemon/unionfind.h>

namespace network_design {
namespace mip {
namespace branch_and_cut {

static char var_name[32];

mip_vars_t::mip_vars_t(GRBModel& model, const instance_t& data)
    : circuit_node(data.graph), circuit_arc(data.graph), star_arc(data.graph) {

    const auto& G = data.graph;
    const int N_PARTITIONS = G.nodeNum();

    // Nós no circuito
    for (Graph::NodeIt i(G); i != lemon::INVALID; ++i) {
        std::snprintf(var_name, sizeof(var_name), "circuit_node[%d]", G.id(i));
        circuit_node[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, var_name);
        circuit_node[i].set(GRB_IntAttr_PoolIgnore, 1);
    }

    // Arcos no circuito
    for (Graph::ArcIt a(G); a != lemon::INVALID; ++a) {
        std::snprintf(var_name, sizeof(var_name), "circuit_arc[%d,%d]",
                      G.id(G.source(a)), G.id(G.target(a)));
        double cost = data.circuit_cost_factor * data.edge_cost[a];

        // Desabilita arestas sem peso definido
        bool usable = cost < 1e99;
        circuit_arc[a] = model.addVar(0.0, usable, cost, GRB_BINARY, var_name);
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

    add_partition_constraints();
    add_star_arc_constraints();
    add_circuit_arc_constraints();

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
                      "Partition %d - Capacity limit", G.id(u));
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
                      "Star node %d In-Degree 1", G.id(v));
        model.addConstr(in_degree_expr == 1 - vars.circuit_node[v],
                        constr_name);
    }
}

void formulation_t::add_circuit_arc_constraints() {
    // Um arco é usado no máximo uma vez, em alguma das categorias de arco.
    for (Graph::EdgeIt e(G); e != lemon::INVALID; ++e) {
        Graph::Arc a = G.arc(G.u(e), G.v(e));
        model.addConstr(vars.star_arc[a] + vars.star_arc[G.oppositeArc(a)] +
                            vars.circuit_arc[a] +
                            vars.circuit_arc[G.oppositeArc(a)] <=
                        1);
    }

    // Uma aresta só pode estar no circuito se ambos seus extremos estão no
    // circuito.
    for (Graph::ArcIt a(G); a != lemon::INVALID; ++a) {
        Graph::Node s = G.source(a), t = G.target(a);
        model.addConstr(vars.circuit_arc[a] <= vars.circuit_node[s]);
        model.addConstr(vars.circuit_arc[a] <= vars.circuit_node[t]);
    }

    // Todo vértice no circuito deve ter grau de entrada e saída 1.
    for (Graph::NodeIt u(G); u != lemon::INVALID; ++u) {
        GRBLinExpr in_degree_expr;
        for (Graph::InArcIt a(G, u); a != lemon::INVALID; ++a) {
            in_degree_expr += vars.circuit_arc[a];
        }
        std::snprintf(constr_name, sizeof(constr_name),
                      "Circuit node %d - In-Degree 1", G.id(u));
        model.addConstr(in_degree_expr == vars.circuit_node[u], constr_name);

        GRBLinExpr out_degree_expr;
        for (Graph::OutArcIt a(G, u); a != lemon::INVALID; ++a) {
            out_degree_expr += vars.circuit_arc[a];
        }
        std::snprintf(constr_name, sizeof(constr_name),
                      "Circuit node %d - Out-Degree 1", G.id(u));
        model.addConstr(out_degree_expr == vars.circuit_node[u], constr_name);
    }
}

void formulation_t::find_violated_integer_cuts() {
    Graph::EdgeMap<double> capacity(G);
    for (Graph::EdgeIt e(G); e != lemon::INVALID; ++e) {
        Graph::Node u = G.u(e), v = G.v(e);
        Graph::Arc a = G.arc(u, v), b = G.arc(v, u);
        capacity[e] = getSolution(vars.circuit_arc[a]) +
                      2 * getSolution(vars.star_arc[a]) +
                      getSolution(vars.circuit_arc[b]) +
                      2 * getSolution(vars.star_arc[b]);
    }

    Graph::NodeMap<int> aux_map(G);
    lemon::UnionFind<Graph::NodeMap<int>> components(aux_map);
    for (Graph::NodeIt v(G); v != lemon::INVALID; ++v) {
        components.insert(v);
    }
    for (Graph::EdgeIt e(G); e != lemon::INVALID; ++e) {
        if (capacity[e] > 1 - 1e-4) {
            components.join(G.u(e), G.v(e));
        }
    }

    int num_components = 0;
    std::vector<int> component_index(G.nodeNum(), -1);
    for (Graph::NodeIt v(G); v != lemon::INVALID; ++v) {
        int component = components.find(v);
        if (component_index[component] < 0) {
            component_index[component] = num_components++;
        }
    }

    for (int i = 0; i < num_components; i++) {
        GRBLinExpr out_expr;
        for (Graph::ArcIt a(G); a != lemon::INVALID; ++a) {
            if (component_index[components.find(G.source(a))] == i &&
                component_index[components.find(G.target(a))] != i) {
                out_expr += vars.circuit_arc[a] + vars.star_arc[a] +
                            vars.star_arc[G.oppositeArc(a)];
            }
        }
        addLazy(out_expr >= 1);
    }
}

void formulation_t::find_violated_fractional_cuts() {
    Graph::EdgeMap<double> capacity(G);
    for (Graph::EdgeIt e(G); e != lemon::INVALID; ++e) {
        Graph::Node u = G.u(e), v = G.v(e);
        Graph::Arc a = G.arc(u, v), b = G.arc(v, u);
        capacity[e] = getSolution(vars.circuit_arc[a]) +
                      2 * getSolution(vars.star_arc[a]) +
                      getSolution(vars.circuit_arc[b]) +
                      2 * getSolution(vars.star_arc[b]);
    }

    std::vector<Graph::Edge> frac_edges, one_edges;
    for (Graph::EdgeIt e(G); e != lemon::INVALID; ++e) {
        if (capacity[e] > 1 - 1e-4) {
            one_edges.push_back(e);
        } else if (capacity[e] > 1e-4) {
            frac_edges.push_back(e);
        }
    }

    Graph::NodeMap<int> aux_map(G);
    lemon::UnionFind<Graph::NodeMap<int>> components(aux_map);
    for (Graph::NodeIt v(G); v != lemon::INVALID; ++v) {
        components.insert(v);
    }
    for (auto e : one_edges) {
        components.join(G.u(e), G.v(e));
    }

    int num_components = 0;
    std::vector<int> component_index(G.nodeNum(), -1);
    for (Graph::NodeIt v(G); v != lemon::INVALID; ++v) {
        int component = components.find(v);
        if (component_index[component] < 0) {
            component_index[component] = num_components++;
        }
    }

    using LGraph = lemon::ListGraph;

    LGraph H;

    std::vector<LGraph::Node> H_index(num_components);
    for (int i = 0; i < num_components; i++) {
        H_index[i] = H.addNode();
    }

    LGraph::EdgeMap<double> H_capacity(H);
    for (auto e : frac_edges) {
        auto u = G.u(e), v = G.v(e);
        auto hu = H_index[component_index[components.find(u)]],
             hv = H_index[component_index[components.find(v)]];
        auto he = H.addEdge(hu, hv);
        H_capacity[he] = capacity[e];
    }

    lemon::GomoryHu<LGraph, LGraph::EdgeMap<double>> gomory_hu(H, H_capacity);
    gomory_hu.run();

    // Controla o número de cortes dependendo da profundidade da busca.
    int n = 0, n_max = (size_t)getDoubleInfo(GRB_CB_MIPSOL_NODCNT) / 2;

    LGraph::NodeMap<bool> cutmap(H, false);
    for (LGraph::NodeIt u(H); u != lemon::INVALID && n <= n_max; ++u) {
        if (gomory_hu.predNode(u) == lemon::INVALID ||
            gomory_hu.predValue(u) > 2.0 - 1e-5) {
            continue;
        }

        gomory_hu.minCutMap(u, gomory_hu.predNode(u), cutmap);

        GRBLinExpr out_expr;
        for (Graph::ArcIt a(G); a != lemon::INVALID; ++a) {
            auto s = G.source(a), t = G.target(a);
            auto hs = H_index[component_index[components.find(s)]],
                 ht = H_index[component_index[components.find(t)]];
            if (cutmap[hs] && !cutmap[ht]) {
                out_expr += vars.circuit_arc[a] + vars.star_arc[a] +
                            vars.star_arc[G.oppositeArc(a)];
            }
        }

        addLazy(out_expr >= 1);
        n++;
    }
}

void formulation_t::find_violated_blossom() {
    Graph::EdgeMap<double> capacity(G);
    for (Graph::EdgeIt e(G); e != lemon::INVALID; ++e) {
        Graph::Node u = G.u(e), v = G.v(e);
        Graph::Arc a = G.arc(u, v), b = G.arc(v, u);
        capacity[e] =
            getNodeRel(vars.circuit_arc[a]) + 2 * getNodeRel(vars.star_arc[a]) +
            getNodeRel(vars.circuit_arc[b]) + 2 * getNodeRel(vars.star_arc[b]);
    }

    lemon::GomoryHu<Graph, Graph::EdgeMap<double>> gomory_hu(G, capacity);
    gomory_hu.run();

    Graph::NodeMap<int> in_degree(G);
    Graph::NodeMap<int> cut_size(G);

    std::list<Graph::Node> sources, topological_sort;
    for (Graph::NodeIt v(G); v != lemon::INVALID; ++v) {
        in_degree[v] = 0;
        cut_size[v] = 1;
    }

    for (Graph::NodeIt u(G); u != lemon::INVALID; ++u) {
        auto v = gomory_hu.predNode(u);
        if (v == lemon::INVALID) {
            continue;
        }
        in_degree[v]++;
    }

    for (Graph::NodeIt v(G); v != lemon::INVALID; ++v) {
        if (in_degree[v] == 0) {
            sources.push_back(v);
        }
    }

    while (!sources.empty()) {
        auto u = sources.front();
        sources.pop_front();
        topological_sort.push_back(u);
        auto v = gomory_hu.predNode(u);
        if (v != lemon::INVALID) {
            in_degree[v]--;
            if (in_degree[v] == 0) {
                sources.push_back(v);
            }
        }
    }

    while (!topological_sort.empty()) {
        auto u = topological_sort.front();
        topological_sort.pop_front();
        auto v = gomory_hu.predNode(u);
        if (v != lemon::INVALID) {
            cut_size[v] = cut_size[v] + cut_size[u];
        }
    }

    for (Graph::NodeIt u(G); u != lemon::INVALID; ++u) {
        GRBLinExpr expr;
        if (cut_size[u] % 2 == 0 || gomory_hu.predValue(u) > 1.0 - 1e-5) {
            continue;
        }

        for (decltype(gomory_hu)::MinCutEdgeIt e(gomory_hu, u,
                                                 gomory_hu.predNode(u));
             e != lemon::INVALID; ++e) {
            Graph::Arc a = G.arc(G.u(e), G.v(e)), b = G.oppositeArc(a);
            expr += vars.circuit_arc[a] + vars.circuit_arc[b] +
                    2 * (vars.star_arc[a] + vars.star_arc[b]);
        }

        addCut(expr >= 1.0);
    }
}

void formulation_t::callback() {
    try {
        switch (where) {
        case GRB_CB_MIPSOL:
            find_violated_fractional_cuts();
            break;

        case GRB_CB_MIPNODE:
            if (getIntInfo(GRB_CB_MIPNODE_STATUS) == GRB_OPTIMAL) {
                // find_violated_blossom();
                // find_violated_fractional_cuts();
            }
            break;

        default:
            break;
        }
    } catch (const GRBException& ex) {
        std::cerr << "ERROR: " << ex.getMessage() << " (" << ex.getErrorCode()
                  << ")" << std::endl;
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
                      << " [color=blue, label=" << data.edge_cost[a] << "];"
                      << std::endl;
            std::cout << "\t" << G.id(t) << " -> " << G.id(s)
                      << " [color=blue, style=dashed];" << std::endl;
        }
    }
    for (Graph::ArcIt a(G); a != lemon::INVALID; ++a) {
        Graph::Node s = G.source(a), t = G.target(a);
        if (vars.circuit_arc[a].get(GRB_DoubleAttr_X) >= 0.5) {
            std::cout << "\t" << G.id(s) << " -> " << G.id(t)
                      << " [color=red, label="
                      << data.circuit_cost_factor * data.edge_cost[a] << "];"
                      << std::endl;
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

} // namespace branch_and_cut
} // namespace mip
} // namespace network_design
