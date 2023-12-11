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
     * `star_edge[e]` = 1 se a aresta `e` está em alguma estrela, 0 caso
     *                  contrário.
     */
    Graph::EdgeMap<GRBVar> star_edge;

    mip_vars_t(GRBModel&, const instance_t&);
};

/**
 * Constrói uma solução a partir das variáveis da formulação polinomial.
 *
 * Deve ser chamado após a função `solve` da formulação correspondente.
 */
solution_t build_solution(const instance_t&, const mip_vars_t&);

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
    void add_circuit_arc_constraints();

    /**
     * Adiciona restrições de ordem para os vértices no circuito. Garante que o
     * circuito de fato forma um ciclo (e não um 2-fator qualquer).
     */
    void add_circuit_order_constraints();

    /**
     * Adiciona restrições de quebra de simetria para os vértices no circuito.
     */
    void add_circuit_symmetry_breaking_constraints();

  public:
    formulation_t(const instance_t&, const GRBEnv& env);
    solution_t solve();
};
} // namespace polynomial
} // namespace mip
} // namespace network_design

#endif // _NETWORK_DESIGN_MIP_POLYNOMIAL_HPP
