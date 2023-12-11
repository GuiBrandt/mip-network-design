#ifndef _NETWORK_DESIGN_MIP_BRANCH_AND_CUT_HPP
#define _NETWORK_DESIGN_MIP_BRANCH_AND_CUT_HPP

#include <cassert>
#include <vector>

#include <gurobi_c++.h>

#include "../defs.hpp"

namespace network_design {
namespace mip {
namespace branch_and_cut {

/**
 * Estrutura para as variáveis da formulação polinomial do problema.
 */
struct mip_vars_t {
    using NodePartitionVars = Graph::NodeMap<std::vector<GRBVar>>;

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
    NodePartitionVars circuit_partition_node;

    /**
     * Variáveis de incidência dos nós do circuito.
     *
     * `circuit_node[i]` = 1 se o nó `i` está no circuito, 0 caso contrário.
     */
    Graph::NodeMap<GRBVar> circuit_node;

    /**
     * Variáveis de incidência das arestas do circuito.
     *
     * `circuit_edge[a]` = 1 se a aresta `a` está no circuito, 0 caso contrário.
     */
    Graph::EdgeMap<GRBVar> circuit_edge;

    /**
     * Variáveis de incidência das arestas das estrelas.
     *
     * `star_edge[e]` = 1 se a aresta `e` está em alguma estrela, 0 caso
     *                  contrário.
     */
    Graph::EdgeMap<GRBVar> star_edge;

    mip_vars_t(GRBModel&, const instance_t&);
};

/**
 * Classe para a formulação polinomial do problema.
 */
class formulation_t : public GRBCallback {
  private:
    const instance_t& instance;
    GRBModel model;
    mip_vars_t vars;

    const Graph& G;
    const int N_PARTITIONS;

    /**
     * Adiciona restrições forçando os valores das variáveis de uso das
     * partições.
     */
    void add_used_partition_constraints();

    /**
     * Adiciona restrições forçando as partições de menor índice a serem usadas
     * para evitar simetrias.
     */
    void add_partition_symmetry_breaking_constraints();

    /**
     * Adiciona restrições de empacotamento dos vértices nas partições.
     */
    void add_partition_packing_constraints();

    /**
     * Adiciona restrições de arestas de estrela.
     *
     * Uma aresta deve estar em uma estrela se ambos seus extremos estão na
     * mesma partição e um deles está no circuito.
     */
    void add_star_edge_constraints();

    /**
     * Adiciona restrições para vértices serem escolhidos no circuito.
     */
    void add_circuit_node_constraints();

    /**
     * Adiciona restrições para arcos no circuito.
     */
    void add_circuit_edge_constraints();

    /**
     * Adiciona restrições de quebra de simetria para os vértices no circuito.
     */
    void add_circuit_symmetry_breaking_constraints();

  protected:
    void callback() override;

  public:
    formulation_t(const instance_t&, const GRBEnv& env);
    solution_t solve();
};

} // namespace branch_and_cut
} // namespace mip
} // namespace network_design

#endif // _NETWORK_DESIGN_MIP_BRANCH_AND_CUT_HPP
