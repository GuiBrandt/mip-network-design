#ifndef _NETWORK_DESIGN_MIP_BASE_HPP
#define _NETWORK_DESIGN_MIP_BASE_HPP

#include <cassert>
#include <vector>

#include <gurobi_c++.h>

#include "../defs.hpp"

namespace network_design {
namespace mip {

/**
 * Estrutura comum para as variáveis das formulações do problema.
 */
struct mip_vars_base {
    /**
     * Variáveis de incidência dos nós do circuito.
     *
     * `circuit_node[i]` = 1 se o nó `i` está no circuito, 0 caso contrário.
     */
    Graph::NodeMap<GRBVar> circuit_node;

    /**
     * Variáveis de incidência dos arcos do circuito.
     *
     * `circuit_arc[a]` = 1 se a aresta `a` está no circuito, 0 caso contrário.
     */
    Graph::ArcMap<GRBVar> circuit_arc;

    /**
     * Variáveis de incidência dos arcos das estrelas.
     *
     * `star_arc[a]` = 1 se o arco `a` está em alguma estrela, 0 caso
     *                  contrário.
     */
    Graph::ArcMap<GRBVar> star_arc;

    mip_vars_base(GRBModel&, const instance_t&);
};

/**
 * Estrutura comum para as formulações do problema.
 */
template <typename Vars> class formulation_base {
    static_assert(std::is_base_of<mip_vars_base, Vars>::value,
                  "variable type must be derived from `mip_vars_base`");

  protected:
    const instance_t& instance;

    GRBModel model;
    Vars vars;

  private:
    /**
     * Adiciona restrições de empacotamento dos vértices.
     */
    void add_partition_constraints() {
        static char constr_name[32];

        const auto& G = instance.graph;

        // Limite de capacidade.
        for (Graph::NodeIt u(G); u != lemon::INVALID; ++u) {
            GRBLinExpr expr;
            for (Graph::OutArcIt a(G, u); a != lemon::INVALID; ++a) {
                expr += instance.node_weight[G.target(a)] * vars.star_arc[a];
            }
            std::snprintf(constr_name, sizeof(constr_name),
                          "partition_capacity[%d]", G.id(u));
            model.addConstr(expr <=
                                (instance.capacity - instance.node_weight[u]) *
                                    vars.circuit_node[u],
                            constr_name);
        }

        // Pelo menos 3 partes são usadas.
        GRBLinExpr expr;
        for (Graph::NodeIt u(G); u != lemon::INVALID; ++u) {
            expr += vars.circuit_node[u];
        }
        model.addConstr(expr >= 3, "At least 3 parts");
    }

    /**
     * Adiciona restrições para arcos no circuito.
     */
    void add_basic_topology_constraints() {
        static char constr_name[32];

        const auto& G = instance.graph;

        // Exatamente um arco de estrela entra num nó que não está no circuito,
        // e nenhum entra num nó que está no circuito.
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
            model.addConstr(in_degree_expr == vars.circuit_node[u],
                            constr_name);

            std::snprintf(constr_name, sizeof(constr_name),
                          "circuit_node_out_degree[%d]", G.id(u));
            model.addConstr(out_degree_expr == vars.circuit_node[u],
                            constr_name);
        }
    }

    solution_t build_solution() {
        const auto& G = instance.graph;

        solution_t solution(instance);

        Graph::NodeMap<int> circuit_partition_index(G);

        // Nós do circuito de acordo com a ordem do circuito
        Graph::Node current, prev;
        for (Graph::NodeIt v(G); v != lemon::INVALID; ++v) {
            if (vars.circuit_node[v].get(GRB_DoubleAttr_X) >= 0.5) {
                current = v;
                break;
            }
        }
        int max_part = -1;
        do {
            circuit_partition_index[current] = ++max_part;
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

        // Partição dos nós
        for (Graph::NodeIt v(G); v != lemon::INVALID; ++v) {
            Graph::Node center;
            if (vars.circuit_node[v].get(GRB_DoubleAttr_X) >= 0.5) {
                center = v;
            } else {
                for (Graph::InArcIt a(G, v); a != lemon::INVALID; ++a) {
                    if (vars.star_arc[a].get(GRB_DoubleAttr_X) >= 0.5) {
                        center = G.source(a);
                        break;
                    }
                }
            }
            solution.partition[v] = circuit_partition_index[center];
        }

        return solution;
    }

  public:
    formulation_base(const instance_t& instance, const GRBEnv& env)
        : instance(instance), model(env), vars(model, instance) {
        model.set(GRB_StringAttr_ModelName, "Um Problema de Network Design");
        model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
        model.update();
        add_partition_constraints();
        add_basic_topology_constraints();
    };

    solution_t solve() {
        model.optimize();
        return build_solution();
    }
};

} // namespace mip
} // namespace network_design

#include "../defs.hpp"

#endif // _NETWORK_DESIGN_MIP_BASE_HPP
