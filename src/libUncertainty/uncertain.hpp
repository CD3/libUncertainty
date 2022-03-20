#pragma once

#include <array>
#include <cmath>
#include <iostream>
#include <type_traits>
#include <utility>

#include "./statistics.hpp"
#include "./utils.hpp"
#include "./tags.hpp"

/** @file uncertain.hpp
 * @brief The uncertain class
 * @author C.D. Clark III
 * @date 03/05/22
 */

namespace libUncertainty
{
template<typename NT, typename UT = decltype(NT() - NT())>
class uncertain
{
 public:
  using nominal_type     = NT;
  using uncertainty_type = UT;
  using this_type        = uncertain<NT, UT>;

  uncertain()                 = default;
  uncertain(const uncertain&) = default;
  uncertain(uncertain&)       = default;
  uncertain(uncertain&&)      = default;
  uncertain& operator=(const uncertain&) = default;
  uncertain& operator=(uncertain&&) = default;

  uncertain(nominal_type nom) : m_storage(std::move(nom), 0) {}
  uncertain(nominal_type nom, uncertainty_type unc) : m_storage(std::move(nom), std::move(unc)) {}
  template<typename T, typename std::enable_if<std::is_same<decltype(std::declval<T>()[0], int()), int>::value, int>::type = 0>
  uncertain(const T& data) : m_storage(data[0], data[1])
  {
  }

  nominal_type     nominal() const { return m_storage.first; }
  uncertainty_type uncertainty() const { return m_storage.second; }

  void nominal(nominal_type a_val) { m_storage.first = a_val; }
  void uncertainty(uncertainty_type a_val) { m_storage.second = a_val; }

  nominal_type upper() const { return m_storage.first + static_cast<nominal_type>(m_storage.second); }
  nominal_type lower() const { return m_storage.first - static_cast<nominal_type>(m_storage.second); }

  auto relative_uncertainty() const -> decltype(this->uncertainty() / this->nominal()) { return this->uncertainty() / this->nominal(); }

  friend std::ostream& operator<<(std::ostream& out, const uncertain<nominal_type, uncertainty_type>& a_val)
  {
    out << a_val.nominal() << " +/- " << a_val.uncertainty();
    return out;
  }

  this_type normalize(size_t n = 1) const
  {
    uncertainty_type unc = sigfig_round(this->uncertainty(), n);
    // the nominal value should be rounded to the same decimal position
    // as the last significant figure of uncertainty.
    // this means that the nominal value and the uncertainty can have
    // a different number of significant figures.
    // it turns out, that the difference will be equal the difference in
    // the exponents of the two numbers expressed in scientific notation.
    // NOTE: We are converting the uncertainty value to the nominal_type here
    // to handle cases when the two are given in different units.
    int          unc_exp = scientific_notation_exponent(nominal_type(unc));
    int          nom_exp = scientific_notation_exponent(this->nominal());
    nominal_type nom     = sigfig_round(this->nominal(), n + nom_exp - unc_exp);
    return {nom, unc};
  }

 private:
  std::pair<nominal_type, uncertainty_type> m_storage;
};

template<typename T>
auto make_uncertain(const T& a_arg)
-> decltype(make_uncertain(a_arg, priority<2>{}))
{
  return make_uncertain(a_arg, priority<2>{});
}

template<typename T1, typename T2>
auto make_uncertain(const T1& a_arg1, const T2& a_arg2)
-> decltype(make_uncertain(a_arg1, a_arg2, priority<2>{}))
{
  return make_uncertain(a_arg1, a_arg2, priority<2>{});
}

template<typename NT>
uncertain<NT> make_uncertain(const NT& nom, priority<0>)
{
  return uncertain<NT>(nom);
}

template<typename NT, typename UT>
uncertain<NT, UT> make_uncertain(const NT& nom, const UT& unc, priority<0>)
{
  return uncertain<NT, UT>(nom, unc);
}

/**
 * Make an uncertain value from a sequance of values using the standard error of the mean for the uncertainty.
 */
template<typename iterator>
auto make_uncertain(iterator begin, iterator end, priority<1>)
-> decltype(make_uncertain(average(begin, end), standard_error_of_the_mean(begin, end), priority<0>{}))
{
  return make_uncertain(average(begin, end), standard_error_of_the_mean(begin, end),priority<0>{});
}

/**
 * Make an uncertain value from a sequence of values. This overload uses the standard deviation
 * for the uncertainty.
 *
 * @param begin Iterator pointing to beginning of sequence.
 * @param end Iterator pointing to end of sequence.
 * @param tag instance of tags::use_stdev_for_error used to select this overload. To use, pass `tags::use_stdev_for_error{}` as the third argument.
 * @param degrees_of_freedom_reduction Reduce the number of degrees of freedom when computing the standard deviation. Default is 1.
 */
template<typename iterator>
auto make_uncertain(iterator begin, iterator end, tags::use_stdev_for_error tag, size_t degrees_of_freedom_reduction = 1)
-> decltype(make_uncertain(average(begin, end), standard_deviation(begin, end, degrees_of_freedom_reduction), priority<0>{}))
{
  return make_uncertain(average(begin, end), standard_deviation(begin, end),priority<0>{});
}

}  // namespace libUncertainty
