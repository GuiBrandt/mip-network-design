#include <stdint.h>

#include <numeric>
#include <random>
#include <unordered_map>
#include <unordered_set>

#include <gurobi_c++.h>
#include <lemon/full_graph.h>

// problem_data read_instance();

using Graph = lemon::FullGraph;

struct problem_data {
    Graph graph;
    Graph::EdgeMap<uint64_t> edge_cost;
    Graph::NodeMap<uint64_t> node_weight;
    uint64_t capacity;
    uint64_t circuit_cost_factor;

    problem_data(const Graph& graph)
        : graph(graph), edge_cost(graph), node_weight(graph) {}
};

using solution = std::vector<Graph::Edge>;

solution solve(const GRBEnv& env, const problem_data& data) {
    GRBModel model(env);
    model.set(GRB_StringAttr_ModelName, "Um Problema de Network Design");
    model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);

    // Número de partições (no pior caso, igual ao número de vértices).
    auto N = data.graph.nodeNum();

    //-------------------------------------------------------------------------
    // Variáveis
    //-------------------------------------------------------------------------
    using NodePartitionVars = Graph::NodeMap<std::vector<GRBVar>>;
    using EdgePartitionVars = Graph::EdgeMap<std::vector<GRBVar>>;

    char var_name[32];

    // Variáveis de incidência das partições.
    // partition[i][j] = 1 se o nó i está na partição j, 0 caso contrário.
    NodePartitionVars partition(data.graph);
    for (Graph::NodeIt i(data.graph); i != lemon::INVALID; ++i) {
        partition[i].resize(N);
        for (int j = 0; j < N; j++) {
            std::snprintf(var_name, sizeof(var_name), "partition[%d,%d]",
                          data.graph.id(i), j);
            partition[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, var_name);
        }
    }

    // Variáveis indicadoras de uso das partições.
    // partition_used[j] = 1 se a partição j tem algum nó, 0 caso contrário.
    std::vector<GRBVar> partition_used(N);
    for (int j = 0; j < N; j++) {
        std::snprintf(var_name, sizeof(var_name), "partition_used[%d]", j);
        partition_used[j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, var_name);
    }

    // Variáveis de incidência dos nós do circuito.
    // circuit_node[i][j] = 1 se o nó i está no circuito e na partição j, 0
    //                      caso contrário.
    NodePartitionVars circuit_node(data.graph);
    for (Graph::NodeIt i(data.graph); i != lemon::INVALID; ++i) {
        circuit_node[i].resize(N);
        for (int j = 0; j < N; j++) {
            std::snprintf(var_name, sizeof(var_name), "circuit_node[%d,%d]",
                          data.graph.id(i), j);
            circuit_node[i][j] =
                model.addVar(0.0, 1.0, 0.0, GRB_BINARY, var_name);
        }
    }

    // Ordem no circuito de um nó.
    Graph::NodeMap<GRBVar> circuit_order(data.graph);
    for (Graph::NodeIt i(data.graph); i != lemon::INVALID; ++i) {
        std::snprintf(var_name, sizeof(var_name), "circuit_order[%d]",
                      data.graph.id(i));
        circuit_order[i] = model.addVar(0.0, N - 1, 0.0, GRB_INTEGER, var_name);
    }

    // Variáveis de incidência dos arcos do circuito.
    // circuit_edge[u][v] = 1 se o arco (u, v) está no circuito, 0 caso
    // contrário.
    Graph::ArcMap<GRBVar> circuit_edge(data.graph);
    for (Graph::ArcIt e(data.graph); e != lemon::INVALID; ++e) {
        Graph::Node u = data.graph.source(e), v = data.graph.target(e);
        std::snprintf(var_name, sizeof(var_name), "circuit_edge[%d,%d]",
                      data.graph.id(u), data.graph.id(v));
        uint64_t cost =
            data.circuit_cost_factor * data.edge_cost[data.graph.edge(u, v)];
        circuit_edge[e] = model.addVar(0.0, 1.0, cost, GRB_BINARY, var_name);
    }

    // Variáveis de incidência das arestas das estrelas.
    // star_edge[e][j] = 1 se a aresta e está na estrela j, 0 caso contrário.
    EdgePartitionVars star_edge(data.graph);
    for (Graph::EdgeIt e(data.graph); e != lemon::INVALID; ++e) {
        star_edge[e].resize(N);
        for (int j = 0; j < N; j++) {
            std::snprintf(var_name, sizeof(var_name), "star_edge[%d,%d]",
                          data.graph.id(e), j);
            uint64_t cost = data.edge_cost[e];
            star_edge[e][j] =
                model.addVar(0.0, 1.0, cost, GRB_BINARY, var_name);
        }
    }

    //-------------------------------------------------------------------------
    // Restrições
    //-------------------------------------------------------------------------
    char constr_name[64];

    // Restrições das variáveis de uso das partições.
    for (int j = 0; j < N; j++) {
        GRBLinExpr expr;
        for (Graph::NodeIt i(data.graph); i != lemon::INVALID; ++i) {
            std::snprintf(constr_name, sizeof(constr_name),
                          "Partition %d - Used with node %d", j,
                          data.graph.id(i));
            model.addConstr(partition_used[j] >= partition[i][j], constr_name);
        }

        // Força as partições de menor índice a serem usadas para evitar
        // simetrias.
        if (j > 0) {
            std::snprintf(constr_name, sizeof(constr_name),
                          "Partition %d - Symmetry breaking", j);
            model.addConstr(partition_used[j] <= partition_used[j - 1],
                            constr_name);
        }
    }

    // Restrições de empacotamento nas partições.
    for (int j = 0; j < N; j++) {
        GRBLinExpr expr;
        for (Graph::NodeIt i(data.graph); i != lemon::INVALID; ++i) {
            expr += partition[i][j];
        }
        std::snprintf(constr_name, sizeof(constr_name),
                      "Partition %d - Capacity limit", j);
        model.addConstr(expr <= data.capacity, constr_name);
    }
    for (Graph::NodeIt i(data.graph); i != lemon::INVALID; ++i) {
        GRBLinExpr expr;
        for (int j = 0; j < N; j++) {
            expr += partition[i][j];
        }
        std::snprintf(constr_name, sizeof(constr_name),
                      "Node %d - Exactly one partition", data.graph.id(i));
        model.addConstr(expr == 1, constr_name);
    }

    // Restrições de vértices no circuito.
    for (int j = 0; j < N; j++) {
        // Exatamente um vértice por partição usada está no circuito.
        GRBLinExpr expr;
        for (Graph::NodeIt i(data.graph); i != lemon::INVALID; ++i) {
            expr += circuit_node[i][j];
        }
        std::snprintf(constr_name, sizeof(constr_name),
                      "Partition %d - Exactly one circuit node", j);
        model.addConstr(expr == partition_used[j], constr_name);

        // O vértice só pode estar no circuito em uma partição se de fato
        // estiver na partição.
        for (Graph::NodeIt i(data.graph); i != lemon::INVALID; ++i) {
            std::snprintf(constr_name, sizeof(constr_name),
                          "Partition %d - Circuit node %d validity", j,
                          data.graph.id(i));
            model.addConstr(circuit_node[i][j] <= partition[i][j], constr_name);
        }
    }

    // O vértice na primeira partição tem ID menor que todos os outros vértices
    // no ciclo. Elimina simetrias de rotação do ciclo.
    for (Graph::EdgeIt e(data.graph); e != lemon::INVALID; ++e) {
        Graph::Node u = data.graph.u(e), v = data.graph.v(e);
        if (data.graph.id(u) > data.graph.id(v)) {
            std::swap(u, v);
        }
        GRBLinExpr is_circuit_node_expr;
        for (int j = 0; j < N; j++) {
            is_circuit_node_expr += circuit_node[u][j];
        }
        model.addConstr(circuit_node[u][0] >=
                            circuit_node[v][0] - (1 - is_circuit_node_expr),
                        constr_name);
    }

    // Restrições de arestas de estrela. Uma aresta deve estar em uma estrela
    // se ambos seus extremos estão na mesma partição e um deles está no
    // circuito.
    for (int j = 0; j < N; j++) {
        for (Graph::EdgeIt e(data.graph); e != lemon::INVALID; ++e) {
            Graph::Node u = data.graph.u(e), v = data.graph.v(e);
            std::snprintf(constr_name, sizeof(constr_name),
                          "Partition %d - Star edge {%d, %d}", j,
                          data.graph.id(u), data.graph.id(v));
            model.addConstr(star_edge[e][j] >=
                                circuit_node[u][j] - 2 * (1 - partition[v][j]),
                            constr_name);
            model.addConstr(star_edge[e][j] >=
                                circuit_node[v][j] - 2 * (1 - partition[u][j]),
                            constr_name);
        }
    }

    for (Graph::NodeIt u(data.graph); u != lemon::INVALID; ++u) {
        for (Graph::NodeIt v(data.graph); v != lemon::INVALID; ++v) {
            if (u != v) {
                model.addConstr(circuit_edge[data.graph.arc(u, v)] +
                                    circuit_edge[data.graph.arc(v, u)] <=
                                1);
            }
        }
    }

    // Todo vértice no circuito deve ter grau de entrada e saída igual a 1.
    for (Graph::NodeIt u(data.graph); u != lemon::INVALID; ++u) {
        GRBLinExpr in_degree_expr;
        GRBLinExpr out_degree_expr;
        for (Graph::NodeIt v(data.graph); v != lemon::INVALID; ++v) {
            if (u != v) {
                in_degree_expr += circuit_edge[data.graph.arc(v, u)];
                out_degree_expr += circuit_edge[data.graph.arc(u, v)];
            }
        }

        GRBLinExpr is_circuit_node_expr;
        for (int j = 0; j < N; j++) {
            is_circuit_node_expr += circuit_node[u][j];
        }
        std::snprintf(constr_name, sizeof(constr_name),
                      "Circuit node %d - In-Degree 1", data.graph.id(u));
        model.addConstr(in_degree_expr == is_circuit_node_expr, constr_name);

        std::snprintf(constr_name, sizeof(constr_name),
                      "Circuit node %d - Out-Degree 1", data.graph.id(u));
        model.addConstr(out_degree_expr == is_circuit_node_expr, constr_name);
    }

    // Limite da ordem cíclica.
    GRBLinExpr expr;
    for (int j = 0; j < N; j++) {
        expr += partition_used[j];
    }
    for (Graph::NodeIt i(data.graph); i != lemon::INVALID; ++i) {
        std::snprintf(constr_name, sizeof(constr_name),
                      "Circuit order of node %d - Number of partitions",
                      data.graph.id(i));
        model.addConstr(circuit_order[i] <= expr - 1);
    }

    // Eliminação de subciclo: se uma aresta está no circuito, então um de seus
    // extremos vem antes do outro na ordem cíclica, exceto para o nó na
    // primeira partição.
    for (Graph::ArcIt a(data.graph); a != lemon::INVALID; ++a) {
        std::snprintf(constr_name, sizeof(constr_name),
                      "Circuit arc %d - Order", data.graph.id(a));
        Graph::Node u = data.graph.source(a), v = data.graph.target(a);
        model.addConstr(circuit_order[v] >= circuit_order[u] + 1 -
                                                N * (1 - circuit_edge[a]) -
                                                N * circuit_node[v][0],
                        constr_name);
    }

    // Uma aresta só pode estar no circuito se ambos seus extremos estão no
    // circuito.
    for (Graph::ArcIt a(data.graph); a != lemon::INVALID; ++a) {
        Graph::Node u = data.graph.source(a), v = data.graph.target(a);
        std::snprintf(constr_name, sizeof(constr_name),
                      "Circuit arc %d - Both ends in circuit",
                      data.graph.id(a));
        GRBLinExpr expr;
        for (int j = 0; j < N; j++) {
            expr += circuit_node[u][j] + circuit_node[v][j];
        }
        model.addConstr(2 * circuit_edge[a] <= expr, constr_name);
    }

    model.update();
    model.optimize();
    model.write("model.lp");
    model.write("model.sol");

    std::cout << "digraph {" << std::endl;

    std::vector<Graph::Edge> solution;

    std::cout << "subgraph Circuit {" << std::endl;
    std::cout << "edge [color=red];" << std::endl;
    for (Graph::ArcIt a(data.graph); a != lemon::INVALID; ++a) {
        if (circuit_edge[a].get(GRB_DoubleAttr_X) >= 0.5) {
            Graph::Node s = data.graph.source(a), t = data.graph.target(a);
            std::cout << data.graph.id(s) << " -> " << data.graph.id(t)
                      << "[penwidth="
                      << data.edge_cost[data.graph.edge(s, t)] / 2.0
                      << "];" << std::endl;
            solution.push_back(data.graph.edge(s, t));
        }
    }
    std::cout << "}" << std::endl;

    for (int j = 0; j < N; j++) {
        std::cout << "subgraph Star" << j << " {" << std::endl;
        std::cout << "edge [color=blue];" << std::endl;
        for (Graph::EdgeIt e(data.graph); e != lemon::INVALID; ++e) {
            if (star_edge[e][j].get(GRB_DoubleAttr_X) >= 0.5) {
                Graph::Node u = data.graph.u(e), v = data.graph.v(e);
                std::cout << data.graph.id(u) << " -> " << data.graph.id(v)
                          << "[dir=none, penwidth=" << data.edge_cost[e] / 2.0
                          << "];" << std::endl;
                solution.push_back(e);
            }
        }
        std::cout << "}" << std::endl;
    }

    std::cout << "}" << std::endl;

    return solution;
}

int main(int argc, char* argv[]) {
    std::random_device rd;
    int seed = 1729;

    GRBEnv env;
    env.set(GRB_IntParam_Seed, seed);

    lemon::FullGraph g(20);

    problem_data instance(g);
    instance.capacity = 15;
    instance.circuit_cost_factor = 3;

    std::minstd_rand rng(seed);
    std::uniform_int_distribution<int> weight_dist(1, 10);
    std::uniform_int_distribution<int> cost_dist(1, 10);

    for (Graph::NodeIt v(g); v != lemon::INVALID; ++v) {
        instance.node_weight[v] = weight_dist(rng);
    }

    for (Graph::EdgeIt e(g); e != lemon::INVALID; ++e) {
        instance.edge_cost[e] = cost_dist(rng);
    }

    auto solution = solve(env, instance);
    std::cout << "solution size: " << solution.size() << std::endl;

    return 0;
}
