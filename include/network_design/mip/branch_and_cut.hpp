#ifndef _NETWORK_DESIGN_MIP_BRANCH_AND_CUT_HPP
#define _NETWORK_DESIGN_MIP_BRANCH_AND_CUT_HPP

#include "base.hpp"

namespace network_design {
namespace mip {
namespace branch_and_cut {

using mip_vars_t = mip_vars_base;

/**
 * Classe para a formulação polinomial do problema.
 */
class formulation_t : public formulation_base<mip_vars_t>, public GRBCallback {
  private:
    /**
     * Adiciona restrições de corte violadas de forma preguiçosa.
     */
    void find_violated_integer_cuts();

    /**
     * Adiciona restrições de corte violadas de forma preguiçosa.
     */
    void find_violated_fractional_cuts();

    /**
     * Adiciona restrições blossom violadas de forma preguiçosa.
     */
    void find_violated_blossom();

  protected:
    void callback() override;

  public:
    formulation_t(const instance_t&, const GRBEnv&);
};

} // namespace branch_and_cut
} // namespace mip
} // namespace network_design

#endif // _NETWORK_DESIGN_MIP_BRANCH_AND_CUT_HPP
