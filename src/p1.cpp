#include <stdint.h>

#include <cassert>
#include <fstream>
#include <numeric>
#include <random>
#include <unordered_map>
#include <vector>

#include <gurobi_c++.h>
#include <lemon/full_graph.h>

using Graph = lemon::FullGraph;

struct problem_data {
    /// Grafo da instância.
    const Graph& graph;

    /// Função de custo das arestas.
    Graph::EdgeMap<uint64_t> edge_cost;

    /// Função de peso dos vértices.
    Graph::NodeMap<uint64_t> node_weight;

    /// Capacidade máxima das partições.
    uint64_t capacity;

    /// Fator multiplicador de custo das arestas no circuito (Gama).
    uint64_t circuit_cost_factor;

    problem_data(const Graph& graph)
        : graph(graph), edge_cost(graph), node_weight(graph) {}

    problem_data(const problem_data& other)
        : graph(other.graph), edge_cost(graph), node_weight(graph),
          capacity(other.capacity),
          circuit_cost_factor(other.circuit_cost_factor) {
        for (Graph::EdgeIt e(graph); e != lemon::INVALID; ++e) {
            edge_cost[e] = other.edge_cost[e];
        }
        for (Graph::NodeIt v(graph); v != lemon::INVALID; ++v) {
            node_weight[v] = other.node_weight[v];
        }
    }
};

struct problem_vars {
    using NodePartitionVars = Graph::NodeMap<std::vector<GRBVar>>;
    using EdgePartitionVars = Graph::EdgeMap<std::vector<GRBVar>>;

    /**
     * Variáveis de incidência das partições.
     *
     * `partition[i][j]` = 1 se o nó `i` está na partição `j`, 0 caso contrário.
     */
    NodePartitionVars partition;

    /**
     * Variáveis indicadoras de uso das partições.
     *
     * `partition_used[j]` = 1 se a partição `j` tem algum nó, 0 caso contrário.
     */
    std::vector<GRBVar> partition_used;

    /**
     * Variáveis de incidência dos nós do circuito.
     *
     * `circuit_node[i][j]` = 1 se o nó `i` está no circuito e na partição `j`,
     *                        0 caso contrário.
     */
    NodePartitionVars circuit_node;

    /**
     * Ordem de um nó no circuito.
     *
     * Assume valores de 0 a N - 1, onde N é o número de partições usadas na
     * solução.
     *
     * Para nós fora do circuito, o valor é irrelevante.
     */
    Graph::NodeMap<GRBVar> circuit_order;

    /**
     * Variáveis de incidência dos arcos (considerando a versão orientada do
     * grafo) do circuito.
     *
     * `circuit_arc[a]` = 1 se o arco `a` está no circuito, 0 caso contrário.
     */
    Graph::ArcMap<GRBVar> circuit_arc;

    /**
     * Variáveis de incidência das arestas das estrelas.
     *
     * `star_edge[e][j]` = 1 se a aresta `e` está na estrela `j`, 0 caso
     *                     contrário.
     */
    EdgePartitionVars star_edge;

    problem_vars(GRBModel& model, const problem_data& data)
        : partition(data.graph), partition_used(data.graph.nodeNum()),
          circuit_node(data.graph), circuit_order(data.graph),
          circuit_arc(data.graph), star_edge(data.graph) {
        static char var_name[32];

        const auto& G = data.graph;
        const int N_PARTITIONS = G.nodeNum();

        // Incidência nas partições
        for (Graph::NodeIt i(G); i != lemon::INVALID; ++i) {
            partition[i].resize(N_PARTITIONS);
            for (int j = 0; j < N_PARTITIONS; j++) {
                std::snprintf(var_name, sizeof(var_name), "partition[%d,%d]",
                              G.id(i), j);
                partition[i][j] =
                    model.addVar(0.0, 1.0, 0.0, GRB_BINARY, var_name);
            }
        }

        // Uso das partições
        for (int j = 0; j < N_PARTITIONS; j++) {
            std::snprintf(var_name, sizeof(var_name), "partition_used[%d]", j);
            partition_used[j] =
                model.addVar(0.0, 1.0, 0.0, GRB_BINARY, var_name);
        }

        // Nós no circuito
        for (Graph::NodeIt i(G); i != lemon::INVALID; ++i) {
            circuit_node[i].resize(N_PARTITIONS);
            for (int j = 0; j < N_PARTITIONS; j++) {
                std::snprintf(var_name, sizeof(var_name), "circuit_node[%d,%d]",
                              G.id(i), j);
                circuit_node[i][j] =
                    model.addVar(0.0, 1.0, 0.0, GRB_BINARY, var_name);
            }
        }

        // Ordem dos nós no circuito
        for (Graph::NodeIt i(G); i != lemon::INVALID; ++i) {
            std::snprintf(var_name, sizeof(var_name), "circuit_order[%d]",
                          G.id(i));
            circuit_order[i] =
                model.addVar(0.0, N_PARTITIONS - 1, 0.0, GRB_INTEGER, var_name);
        }

        // Arcos no circuito
        for (Graph::ArcIt e(G); e != lemon::INVALID; ++e) {
            Graph::Node u = G.source(e), v = G.target(e);
            std::snprintf(var_name, sizeof(var_name), "circuit_arc[%d,%d]",
                          G.id(u), G.id(v));
            uint64_t cost =
                data.circuit_cost_factor * data.edge_cost[G.edge(u, v)];
            circuit_arc[e] = model.addVar(0.0, 1.0, cost, GRB_BINARY, var_name);
        }

        // Arestas nas estrelas
        for (Graph::EdgeIt e(G); e != lemon::INVALID; ++e) {
            star_edge[e].resize(N_PARTITIONS);
            for (int j = 0; j < N_PARTITIONS; j++) {
                std::snprintf(var_name, sizeof(var_name), "star_edge[%d,%d]",
                              G.id(e), j);
                uint64_t cost = data.edge_cost[e];
                star_edge[e][j] =
                    model.addVar(0.0, 1.0, cost, GRB_BINARY, var_name);
            }
        }
    }
};

static char constr_name[64];

class problem_constraints_generator {
  private:
    GRBModel& model;
    const problem_vars& vars;
    const problem_data& data;

    const Graph& G;
    const int N_PARTITIONS;

  public:
    problem_constraints_generator(GRBModel& model, const problem_vars& vars,
                                  const problem_data& data)
        : model(model), vars(vars), data(data), G(data.graph),
          N_PARTITIONS(G.nodeNum()) {}

    /**
     * Adiciona restrições forçando os valores das variáveis de uso das
     * partições.
     */
    problem_constraints_generator& add_used_partition_constraints() {
        for (int j = 0; j < N_PARTITIONS; j++) {
            GRBLinExpr expr;
            for (Graph::NodeIt i(G); i != lemon::INVALID; ++i) {
                std::snprintf(constr_name, sizeof(constr_name),
                              "Partition %d - Used with node %d", j, G.id(i));
                model.addConstr(vars.partition_used[j] >= vars.partition[i][j],
                                constr_name);
            }
        }
        return *this;
    }

    /**
     * Adiciona restrições forçando as partições de menor índice a serem usadas
     * para evitar simetrias.
     */
    problem_constraints_generator&
    add_partition_symmetry_breaking_constraints() {
        for (int j = 1; j < N_PARTITIONS; j++) {
            std::snprintf(constr_name, sizeof(constr_name),
                          "Partition %d - Symmetry breaking", j);
            model.addConstr(vars.partition_used[j] <=
                                vars.partition_used[j - 1],
                            constr_name);
        }
        return *this;
    }

    /**
     * Adiciona restrições de empacotamento dos vértices nas partições.
     */
    problem_constraints_generator& add_partition_packing_constraints() {
        // Limite de capacidade.
        for (int j = 0; j < N_PARTITIONS; j++) {
            GRBLinExpr expr;
            for (Graph::NodeIt i(G); i != lemon::INVALID; ++i) {
                expr += data.node_weight[i] * vars.partition[i][j];
            }
            std::snprintf(constr_name, sizeof(constr_name),
                          "Partition %d - Capacity limit", j);
            model.addConstr(expr <= data.capacity, constr_name);
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
        return *this;
    }

    /**
     * Adiciona restrições para vértices serem escolhidos no circuito.
     */
    problem_constraints_generator& add_circuit_node_constraints() {
        for (int j = 0; j < N_PARTITIONS; j++) {
            // Exatamente um vértice da partição está no circuito se ela for
            // usada, ou nenhum, caso contrário.
            GRBLinExpr nodes_in_partition_expr;
            for (Graph::NodeIt i(G); i != lemon::INVALID; ++i) {
                nodes_in_partition_expr += vars.circuit_node[i][j];
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
                model.addConstr(vars.circuit_node[i][j] <= vars.partition[i][j],
                                constr_name);
            }
        }
        return *this;
    }

    /**
     * Adiciona restrições de arestas de estrela.
     *
     * Uma aresta deve estar em uma estrela se ambos seus extremos estão na
     * mesma partição e um deles está no circuito.
     */
    problem_constraints_generator& add_star_edge_constraints() {
        for (int j = 0; j < N_PARTITIONS; j++) {
            for (Graph::EdgeIt e(G); e != lemon::INVALID; ++e) {
                Graph::Node u = G.u(e), v = G.v(e);
                std::snprintf(constr_name, sizeof(constr_name),
                              "Partition %d - Star edge {%d, %d}", j, G.id(u),
                              G.id(v));
                model.addConstr(vars.star_edge[e][j] >=
                                    vars.circuit_node[u][j] -
                                        2 * (1 - vars.partition[v][j]),
                                constr_name);
                model.addConstr(vars.star_edge[e][j] >=
                                    vars.circuit_node[v][j] -
                                        2 * (1 - vars.partition[u][j]),
                                constr_name);
            }
        }
        return *this;
    }

    /**
     * Adiciona restrições para arcos no circuito.
     */
    problem_constraints_generator& add_circuit_arc_constraints() {
        // Um arco só pode estar no circuito se ambos seus extremos estão no
        // circuito.
        for (Graph::ArcIt a(G); a != lemon::INVALID; ++a) {
            Graph::Node s = G.source(a), t = G.target(a);
            std::snprintf(constr_name, sizeof(constr_name),
                          "Circuit arc %d - Both ends in circuit", G.id(a));
            GRBLinExpr expr;
            for (int j = 0; j < N_PARTITIONS; j++) {
                expr += vars.circuit_node[s][j] + vars.circuit_node[t][j];
            }
            model.addConstr(2 * vars.circuit_arc[a] <= expr, constr_name);
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

            GRBLinExpr is_circuit_node_expr;
            for (int j = 0; j < N_PARTITIONS; j++) {
                is_circuit_node_expr += vars.circuit_node[u][j];
            }

            std::snprintf(constr_name, sizeof(constr_name),
                          "Circuit node %d - In-Degree 1", G.id(u));
            model.addConstr(in_degree_expr == is_circuit_node_expr,
                            constr_name);

            std::snprintf(constr_name, sizeof(constr_name),
                          "Circuit node %d - Out-Degree 1", G.id(u));
            model.addConstr(out_degree_expr == is_circuit_node_expr,
                            constr_name);
        }

        return *this;
    }

    /**
     * Adiciona restrições de ordem para os vértices no circuito. Garante que o
     * circuito de fato forma um ciclo (e não um 2-fator qualquer).
     */
    problem_constraints_generator& add_circuit_order_constraints() {
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
                                    N_PARTITIONS * vars.circuit_node[v][0],
                            constr_name);
        }

        return *this;
    }

    /**
     * Adiciona restrições de quebra de simetria para os vértices no circuito.
     */
    problem_constraints_generator& add_circuit_symmetry_breaking_constraints() {
        // O vértice na primeira partição tem ID menor que todos os outros
        // vértices no ciclo. Elimina simetrias de rotação.
        for (Graph::EdgeIt e(G); e != lemon::INVALID; ++e) {
            Graph::Node u = G.u(e), v = G.v(e);

            if (G.id(u) > G.id(v)) {
                std::swap(u, v);
            }

            GRBLinExpr is_circuit_node_expr;
            for (int j = 0; j < N_PARTITIONS; j++) {
                is_circuit_node_expr += vars.circuit_node[u][j];
            }

            model.addConstr(vars.circuit_node[u][0] >=
                                vars.circuit_node[v][0] -
                                    (1 - is_circuit_node_expr),
                            constr_name);
        }
        return *this;
    }
};

struct solution_t {
  private:
    std::vector<Graph::Node> partition_repr;

  public:
    const problem_data& instance;
    const problem_vars& vars;

    /// Partição dos nós do grafo.
    Graph::NodeMap<int> partition;

    /// Nós do circuito, em ordem.
    std::vector<Graph::Node> circuit_nodes;

    solution_t(const problem_data& data, const problem_vars& vars)
        : instance(data), vars(vars), partition(data.graph) {
        const auto& G = data.graph;

        // Partições
        int max_partition = 0;
        for (Graph::NodeIt i(G); i != lemon::INVALID; ++i) {
            for (int j = 0; j < G.nodeNum(); j++) {
                if (vars.partition[i][j].get(GRB_DoubleAttr_X) >= 0.5) {
                    partition[i] = j;
                    if (j > max_partition) {
                        max_partition = j;
                    }
                }
            }
        }

        // Nós do circuito
        circuit_nodes.reserve(max_partition + 1);
        for (Graph::NodeIt v(G); v != lemon::INVALID; ++v) {
            for (int j = 0; j <= max_partition; j++) {
                if (vars.circuit_node[v][j].get(GRB_DoubleAttr_X) >= 0.5) {
                    circuit_nodes.push_back(v);
                    break;
                }
            }
            if (circuit_nodes.size() > max_partition) {
                break;
            }
        }

        assert(circuit_nodes.size() == max_partition + 1);

        // Ordena os nós de acordo com a ordem do circuito
        std::sort(circuit_nodes.begin(), circuit_nodes.end(),
                  [&vars](const Graph::Node& u, const Graph::Node& v) {
                      return vars.circuit_order[u].get(GRB_DoubleAttr_X) <=
                             vars.circuit_order[v].get(GRB_DoubleAttr_X);
                  });

        // Índice reverso para a ordem das partições no circuito
        partition_repr.resize(max_partition + 1);
        for (auto v : circuit_nodes) {
            partition_repr[partition[v]] = v;
        }
    }

    solution_t(const solution_t& other)
        : instance(other.instance), vars(other.vars),
          partition_repr(other.partition_repr), partition(instance.graph),
          circuit_nodes(other.circuit_nodes) {
        for (Graph::NodeIt i(instance.graph); i != lemon::INVALID; ++i) {
            partition[i] = other.partition[i];
        }
    }

    solution_t(const solution_t&& other)
        : instance(other.instance), vars(other.vars),
          partition_repr(std::move(other.partition_repr)),
          partition(instance.graph),
          circuit_nodes(std::move(other.circuit_nodes)) {
        for (Graph::NodeIt i(instance.graph); i != lemon::INVALID; ++i) {
            partition[i] = other.partition[i];
        }
    }

    /// Gera as arestas do circuito na solução.
    std::vector<Graph::Edge> circuit_edges() const {
        const size_t N = circuit_nodes.size();
        std::vector<Graph::Edge> result;
        result.reserve(N);
        for (int i = 0; i < N; i++) {
            result.push_back(instance.graph.edge(circuit_nodes[i],
                                                 circuit_nodes[(i + 1) % N]));
        }
        return result;
    }

    /// Gera as arestas das estrelas na solução.
    std::vector<Graph::Edge> star_edges() const {
        std::vector<Graph::Edge> result;
        result.reserve(instance.graph.nodeNum() - circuit_nodes.size());
        for (Graph::NodeIt u(instance.graph); u != lemon::INVALID; ++u) {
            const auto& v = partition_repr[partition[u]];
            if (u != v) {
                result.push_back(instance.graph.edge(u, v));
            }
        }
        return result;
    }
};

solution_t solve(const GRBEnv& env, const problem_data& data) {
    GRBModel model(env);
    model.set(GRB_StringAttr_ModelName, "Um Problema de Network Design");
    model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);

    for (Graph::EdgeIt e(data.graph); e != lemon::INVALID; ++e) {
        std::cout << data.edge_cost[e] << std::endl;
    }
    for (Graph::NodeIt v(data.graph); v != lemon::INVALID; ++v) {
        std::cout << data.node_weight[v] << std::endl;
    }

    problem_vars vars(model, data);
    problem_constraints_generator(model, vars, data)
        .add_used_partition_constraints()
        .add_partition_packing_constraints()
        .add_partition_symmetry_breaking_constraints()
        .add_star_edge_constraints()
        .add_circuit_node_constraints()
        .add_circuit_arc_constraints()
        .add_circuit_order_constraints()
        .add_circuit_symmetry_breaking_constraints();

    model.update();
    model.optimize();
    model.write("model.lp");
    model.write("model.sol");

    return solution_t(data, vars);
}

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

problem_data random_instance(const Graph& G, int seed) {
    std::minstd_rand rng(seed);

    problem_data instance(G);

    std::uniform_int_distribution<uint64_t> capacity_dist(20, 40);
    instance.capacity = capacity_dist(rng);

    std::uniform_int_distribution<uint64_t> factor_dist(2, 5);
    instance.circuit_cost_factor = factor_dist(rng);

    std::uniform_int_distribution<uint64_t> weight_dist(1, 10);
    for (Graph::NodeIt v(G); v != lemon::INVALID; ++v) {
        instance.node_weight[v] = weight_dist(rng);
    }

    std::uniform_int_distribution<uint64_t> cost_dist(1, 10);
    for (Graph::EdgeIt e(G); e != lemon::INVALID; ++e) {
        instance.edge_cost[e] = cost_dist(rng);
    }

    return instance;
}

int main(int argc, char* argv[]) {
    std::random_device rd;
    int seed = rd();
    std::cout << "seed: " << seed << std::endl;

    std::minstd_rand rng(seed);
    std::uniform_int_distribution<int> grb_seed_dist(0, GRB_MAXINT);

    GRBEnv env;
    env.set(GRB_IntParam_Seed, grb_seed_dist(rng));

    Graph G(16);
    auto instance = random_instance(G, seed);
    auto solution = solve(env, instance);

    view(solution);

    return 0;
}
