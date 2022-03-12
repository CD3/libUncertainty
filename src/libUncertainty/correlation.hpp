#pragma once

/** @file correlation.hpp
  * @brief 
  * @author C.D. Clark III
  * @date 03/11/22
  */

#include <exception>
#include <iostream>
#include <map>

#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "./utils.hpp"

/**
 * A container for storing correlation coefficients in a matrix lay0ut
 */
template<typename T>
struct correlation_matrix {
 private:
  using matrix_type = boost::numeric::ublas::triangular_matrix<T, boost::numeric::ublas::upper>;
  matrix_type m_elements;

 public:
  correlation_matrix() = delete;
  correlation_matrix(size_t N) : m_elements(N, N) {}
  T  operator()(int i, int j) const { return m_elements(std::min(i, j), std::max(i, j)); }
  T& operator()(int i, int j) { return m_elements(std::min(i, j), std::max(i, j)); }

  friend std::ostream& operator<<(std::ostream& out, const correlation_matrix<T>& a_mat)
  {
    out << a_mat.m_elements;
    return out;
  }
};

/**
 * A container for storing a variable with it's correlation coefficients.
 */
namespace libUncertainty
{
template<typename BASE, typename CoeffType = double>
struct add_correlation_coefficient_array : public BASE {
 public:
  using coefficient_type = CoeffType;
  template<typename T>
  using vector_type = std::vector<T>;

  coefficient_type&              get_correlation_coefficient(int i) { return m_correlation_coefficients[i]; }
  coefficient_type               get_correlation_coefficient(int i) const { return m_correlation_coefficients[i]; }
  vector_type<coefficient_type>& get_correlation_coefficients() { return m_correlation_coefficients; }
  vector_type<coefficient_type>  get_correlation_coefficients() const { return m_correlation_coefficients; }

  void add_correlation_coefficient(const coefficient_type& a_coeff) { m_correlation_coefficients.push_back(a_coeff); }
  void set_correlation_coefficient_array_size(size_t N) { m_correlation_coefficients.resize(N); }

  using BASE::BASE;
  using BASE::operator=;

 private:
  vector_type<coefficient_type> m_correlation_coefficients;
};

/**
 * A container for storing correlation coefficients
 */
template<typename T>
struct correlation_store {
 public:
  using key_type = decltype(get_uniq_id());
  using map_type = std::map<key_type, T>;

  void add(const key_type& a_key, const T& a_val)
  {
    if(m_correlation_coefficients.find(a_key) != m_correlation_coefficients.end()) {
      throw std::runtime_error("Correlation entry for " + std::to_string(a_key) + " already exists. Use set(k,v) instead.");
    }
    m_correlation_coefficients[a_key] = a_val;
  }
  void set(const key_type& a_key, const T& a_val) { m_correlation_coefficients[a_key] = a_val; }
  T    get(const key_type& a_key) const { return m_correlation_coefficients.at(a_key); }

 private:
  map_type m_correlation_coefficients;
};
}  // namespace libUncertainty
