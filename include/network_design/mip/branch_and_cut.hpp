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
     * Adiciona restrições de empacotamento dos vértices nas partições.
     */
    void add_partition_constraints();

    /**
     * Adiciona restrições de arestas de estrela.
     *
     * Uma aresta deve estar em uma estrela se ambos seus extremos estão na
     * mesma partição e um deles está no circuito.
     */
    void add_star_arc_constraints();

    /**
     * Adiciona restrições para arcos no circuito.
     */
    void add_circuit_arc_constraints();

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
