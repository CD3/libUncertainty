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
   *
   * DOES NOT HANDLE CORRELATED INPUTS
   */
  template<typename F, typename N1, typename U1, typename N2, typename U2>
  static auto propagate_error(F f, const uncertain<N1, U1>& a1, const uncertain<N2,U2>& a2) -> uncertain<decltype(f(a1.nominal(),a2.nominal()))>
  {
    auto nom = f(a1.nominal(),a2.nominal());
    auto d1 = f(a1.upper(),a2.nominal()) - nom;
    auto d2 = f(a1.nominal(),a2.upper()) - nom;
    auto unc = sqrt(d1*d1 + d2*d2);

    return {nom,unc};
  }

  /**
   * Propagate error through a trinary function f.
   *
   * DOES NOT HANDLE CORRELATED INPUTS
   */
  template<typename F,
    typename N1, typename U1,
    typename N2, typename U2,
    typename N3, typename U3
      >
  static auto propagate_error(F f,
      const uncertain<N1,U1>& a_a1,
      const uncertain<N2,U2>& a_a2,
      const uncertain<N3,U3>& a_a3
      ) -> uncertain<decltype(f(a_a1.nominal(),a_a2.nominal(),a_a3.nominal()))>
  {
    auto nom = f(a_a1.nominal(),a_a2.nominal(),a_a3.nominal());
    auto d1 = f(a_a1.upper(),a_a2.nominal(),a_a3.nominal()) - nom;
    auto d2 = f(a_a1.nominal(),a_a2.upper(),a_a3.nominal()) - nom;
    auto d3 = f(a_a1.nominal(),a_a2.nominal(),a_a3.upper()) - nom;
    auto unc = sqrt(d1*d1 + d2*d2 + d3*d3);

    return {nom,unc};
  }


};
}  // namespace error_propagators

}  // namespace libUncertainty
