#ifndef _NETWORK_DESIGN_MIP_POLYNOMIAL_HPP
#define _NETWORK_DESIGN_MIP_POLYNOMIAL_HPP

#include <cassert>
#include <vector>

#include <gurobi_c++.h>

#include "../defs.hpp"

namespace network_design {
namespace mip {
namespace polynomial {

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
     * Ordem de um nó no circuito.
     *
     * Assume valores de 0 a N - 1, onde N é o número de partições usadas na
     * solução.
     *
     * Para nós fora do circuito, o valor é igual ao nó do circuito em cuja
     * partição o nó está.
     */
    Graph::NodeMap<GRBVar> node_order;

    /**
     * Variáveis de incidência dos arcos (considerando a versão orientada do
     * grafo) do circuito.
     *
     * `circuit_arc[a]` = 1 se o arco `a` está no circuito, 0 caso contrário.
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
class formulation_t {
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

    /**
     * Adiciona restrições de ordem para os vértices. Garante que o grafo é
     * conexo.
     */
    void add_node_order_constraints();

  public:
    formulation_t(const instance_t&, const GRBEnv& env);
    solution_t solve();
};

} // namespace polynomial
} // namespace mip
} // namespace network_design

#endif // _NETWORK_DESIGN_MIP_POLYNOMIAL_HPP
