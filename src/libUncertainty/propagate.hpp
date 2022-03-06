#pragma once
#include "./uncertain.hpp"

/** @file propagate.hpp
  * @brief 
  * @author C.D. Clark III
  * @date 03/05/22
  */

namespace libUncertainty
{
namespace error_propagators
{
struct basic {
  /**
   * Propagate error through a unairy function f.
   */
  template<typename F, typename N1, typename U1>
  static auto propagate_error(F f, const uncertain<N1, U1>& a1) -> uncertain<decltype(f(a1.nominal()))>
  {
    auto nom = f(a1.nominal());
    auto d = f(a1.upper())-nom;
    auto unc = sqrt(d*d);

    return {nom,unc};
  }

  /**
   * Propagate error through a binary function f.
   */
  template<typename F, typename N1, typename U1, typename N2, typename U2>
  static auto propagate_error(F f, const uncertain<N1, U1>& a1, const uncertain<N2,U2>& a2) -> uncertain<decltype(f(a1.nominal(),a2.nominal()))>
  {
    auto nom = f(a1.nominal());
    auto d1 = f(a1.upper()) - nom;
    auto d2 = f(a2.upper()) - nom;
    auto unc = sqrt(d1*d1 + d2*d2);

    return {nom,unc};
  }


};
}  // namespace error_propagators

}  // namespace libUncertainty
