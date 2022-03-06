#pragma once

#include <sstream>
#include <iostream>

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
T sigfig_round(T t, size_t n)
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

}  // namespace libUncertainty
