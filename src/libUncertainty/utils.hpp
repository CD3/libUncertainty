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

/**
 * A function that returns a uniq ID each time it is called.
 */
inline
size_t get_uniq_id()
{
  static size_t id = 0;
  return ++id;
}

/**
 * A mixin class for adding an ID to a variable.
 * Used to track correlations between variables.
 */
template<typename BASE>
class add_id : public BASE
{
  private:
    size_t id = get_uniq_id();

  public:
    using BASE::BASE;
    using BASE::operator=;

    size_t get_id() const {return id;}
    void new_id() {id = get_uniq_id();}
    void clear_id() {id = 0;}


};

template<unsigned P> struct priority: priority<P-1> {};
template<> struct priority<0> {};

template<typename T>
size_t get_id(const T& a_var)
{
  return get_id(a_var,priority<2>{});
}

template<typename T>
size_t get_id(const T& a_var, priority<0>)
{
  return 0;
}

template<typename T>
auto get_id(const T& a_var, priority<1>) -> decltype(a_var.get_id())
{
  return a_var.get_id();
}


template<typename T>
bool is_uncertain(const T& a_var)
{
  return is_uncertain<T>(priority<2>{});
}

template<typename T>
auto is_uncertain(priority<0>) -> decltype(false)
{
  return false;
}

template<typename T>
auto is_uncertain(priority<1>) -> decltype(std::declval<T>().uncertainty(),true)
{
  return true;
}



template<typename T>
auto get_nominal(const T& a_var) -> decltype(get_nominal(a_var,priority<2>{}))
{
  return get_nominal(a_var, priority<2>{});
}

template<typename T>
auto get_nominal(const T& a_var, priority<0>) -> decltype(a_var)
{
  return a_var;
}

template<typename T>
auto get_nominal(const T& a_var, priority<1>) -> decltype(a_var.nominal())
{
  return a_var.nominal();
}




template<typename T>
auto get_uncertainty(const T& a_var) -> decltype(get_uncertainty(a_var,priority<2>{}))
{
  return get_uncertainty(a_var, priority<2>{});
}

template<typename T>
auto get_uncertainty(const T& a_var, priority<0>) -> decltype(T(0))
{
  return T(0);
}

template<typename T>
auto get_uncertainty(const T& a_var, priority<1>) -> decltype(a_var.uncertainty())
{
  return a_var.uncertainty();
}


template<typename T>
auto get_upper(const T& a_var) -> decltype(get_upper(a_var,priority<2>{}))
{
  return get_upper(a_var, priority<2>{});
}

template<typename T>
auto get_upper(const T& a_var, priority<0>) -> decltype(a_var)
{
  return a_var;
}

template<typename T>
auto get_upper(const T& a_var, priority<1>) -> decltype(a_var.upper())
{
  return a_var.upper();
}


template<typename T>
auto get_lower(const T& a_var) -> decltype(get_lower(a_var,priority<2>{}))
{
  return get_lower(a_var, priority<2>{});
}

template<typename T>
auto get_lower(const T& a_var, priority<0>) -> decltype(a_var)
{
  return a_var;
}

template<typename T>
auto get_lower(const T& a_var, priority<1>) -> decltype(a_var.lower())
{
  return a_var.lower();
}

/**
 * A function for getting "zero" for a given type.
 */
template<typename T>
T zero()
{
  return zero<T>(priority<2>{});
}

template<typename T>
T zero(priority<0>)
{
  return static_cast<T>(0);
}

template<typename T>
auto zero(priority<1>) -> decltype( T::from_value(0) )
{
  return T::from_value(0);
}


}  // namespace libUncertainty
