#ifndef _NETWORK_DESIGN_MIP_POLYNOMIAL_HPP
#define _NETWORK_DESIGN_MIP_POLYNOMIAL_HPP

#include "base.hpp"

namespace network_design {
namespace mip {
namespace polynomial {

/**
 * Estrutura para as variáveis da formulação polinomial do problema.
 */
struct mip_vars_t : public mip_vars_base {
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

    mip_vars_t(GRBModel&, const instance_t&);
};

/**
 * Classe para a formulação polinomial do problema.
 */
class formulation_t : public formulation_base<mip_vars_t> {
  private:
    /**
     * Adiciona restrições de ordem para os vértices no circuito. Garante que o
     * grafo é conexo.
     */
    void add_circuit_order_constraints();

  public:
    formulation_t(const instance_t&, const GRBEnv&);
};

} // namespace polynomial
} // namespace mip
} // namespace network_design

#endif // _NETWORK_DESIGN_MIP_POLYNOMIAL_HPP
