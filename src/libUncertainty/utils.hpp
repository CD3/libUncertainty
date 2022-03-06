#pragma once

#include <sstream>
#include <iostream>
#include <cmath>

/** @file utils.hpp
  * @brief Utility functions/classes
  * @author C.D. Clark III
  * @date 03/05/22
  */

namespace libUncertainty
{
  /**
   * Round a number to the given number of significant figures.
   */
template<typename T>
auto sigfig_round(T t, size_t n) -> decltype(std::cin>>t,T())
{
  /**
   * This is a very simple (but probably slow) method.
   * We just use the output stream operator to write
   * the number in scientific notation with n-1 degits after
   * the decimal point, then read it back.
   */
  std::stringstream ss;
  ss.precision(n-1);
  ss << std::scientific << t;
  ss >> t;
  return t;
}

/**
 * Round a Boost.Units quantity<...>-like object to the given number of significant figures.
 */
template<typename T>
auto sigfig_round(T t, size_t n) -> decltype(T::from_value(t.value()),T())
{
  return T::from_value(sigfig_round(t.value(),n));
}

/**
 * Returns the expoenent of a number when expressed in scientific notation.
 *
 * example:
 *
 * 83.521 = 8.3621 x 10^1 ==> return 1
 */
template<typename T>
auto scientific_notation_exponent(T t) -> decltype(std::floor(std::log10(t)),int())
{
  return std::floor(std::log10(t));
}
template<typename T>
auto scientific_notation_exponent(T t) -> decltype(std::floor(std::log10(t.value())),int())
{
  return std::floor(std::log10(t.value()));
}

}  // namespace libUncertainty
