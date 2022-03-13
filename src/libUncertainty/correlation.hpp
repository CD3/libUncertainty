#pragma once

/** @file correlation.hpp
  * @brief 
  * @author C.D. Clark III
  * @date 03/11/22
  */

#include <algorithm>
#include <exception>
#include <iostream>
#include <map>

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/triangular.hpp>

#include "./utils.hpp"

namespace libUncertainty
{
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
  using id_type  = decltype(get_uniq_id());
  using key_type = std::pair<id_type, id_type>;
  using map_type = std::map<key_type, T>;

  key_type make_key(id_type a_id1, id_type a_id2) const
  {
    if(a_id1 > a_id2)
      std::swap(a_id1, a_id2);
    return {a_id1, a_id2};
  }

  /**
   * Add an entry to the correlation store for the id pair (a_id1,a_id2).
   *
   * If an entry already exists, this function will throw a std::runtime_error(...). Use set_with_ids instead.
   */
  void add_with_ids(const id_type& a_id1, const id_type& a_id2, const T& a_val)
  {
    auto key = make_key(a_id1, a_id2);
    if(m_correlation_coefficients.find(key) != m_correlation_coefficients.end()) {
      throw std::runtime_error("Correlation entry for (" + std::to_string(a_id1) + "," + std::to_string(a_id2) + ") already exists. Use set(k,v) instead.");
    }
    m_correlation_coefficients[key] = a_val;
  }

  /**
   * Set an entry in the correlation store for the id pair (a_id1,a_id2).
   *
   * If no entry exists, it is created.
   */
  void set_with_ids(const id_type& a_id1, const id_type& a_id2, const T& a_val)
  {
    m_correlation_coefficients[make_key(a_id1, a_id2)] = a_val;
  }

  /**
   * Get the correlation for the id pari (a_id1,a_id2).
   *
   * If no entry exists, this function returns zero.
   */
  T get_with_ids(const id_type& a_id1, const id_type& a_id2) const
  {
    auto key = make_key(a_id1, a_id2);
    if(m_correlation_coefficients.find(key) == m_correlation_coefficients.end()) {
      return static_cast<T>(0);
    }
    return m_correlation_coefficients.at(key);
  }

  /**
   * Add an entry to the correlation store for a pair of variables.
   *
   * Each variable should have a .get_id() method.
   */
  template<typename U, typename V>
  void add(const U& a_v1, const V& a_v2, const T& a_val)
  {
    add_with_ids(get_id(a_v1), get_id(a_v2), a_val);
  }

  /**
   * Set an entry in the correlation store for a pair of variables.
   *
   * Each variable should have a .get_id() method.
   */
  template<typename U, typename V>
  void set(const U& a_v1, const V& a_v2, const T& a_val)
  {
    set_with_ids(get_id(a_v1), get_id(a_v2), a_val);
  }

  /**
   * Get an entry from the correlation store for a pair of variables.
   *
   * Each variable should have a .get_id() method.
   */
  template<typename U, typename V>
  T get(const U& a_v1, const V& a_v2)
  {
    return get_with_ids(get_id(a_v1), get_id(a_v2));
  }

 private:
  map_type m_correlation_coefficients;
};

/**
 * Returns a reference to a static, global correlation store.
 */
correlation_store<double>& get_global_correlation_store()
{
  static correlation_store<double> store;
  return store;
}
}  // namespace libUncertainty
