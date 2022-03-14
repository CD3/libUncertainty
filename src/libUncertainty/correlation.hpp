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

#include "./utils.hpp"

namespace libUncertainty
{
/**
 * A container for storing correlation coefficients in a matrix lay0ut
 */
template<typename T>
struct correlation_matrix {
 private:
  std::vector<T> m_elements;

  // compute the number of elements that must be stored for an NxN matrix.
  size_t compute_storage_size(size_t a_N) const { return 1 + a_N * (a_N - 1) / 2; }
  // compute the size of a matrix (i.e. N in NxN) from the number of elements stored n.
  size_t compute_matrix_size(size_t a_n) const { return (1 + sqrt(1 + 8 * (a_n - 1))) / 2; }
  // compute the memory index from matrix index
  size_t compute_index(size_t a_i, size_t a_j) const
  {
    // all diagonal elements should be one, and are stored
    // in the last element of the vector
    if(a_i == a_j)
      return m_elements.size() - 1;
    return a_i + a_j - 1;
  }

 public:
  correlation_matrix() = default;
  correlation_matrix(size_t a_N) : m_elements(compute_storage_size(a_N))
  {
    std::fill(m_elements.begin(), m_elements.end() - 1, static_cast<T>(0));
    m_elements[m_elements.size() - 1] = static_cast<T>(1);
  }
  T operator()(int a_i, int a_j) const
  {
    return m_elements[compute_index(a_i, a_j)];
  }
  T& operator()(int a_i, int a_j)
  {
    return m_elements[compute_index(a_i, a_j)];
  }

  friend std::ostream& operator<<(std::ostream& out, const correlation_matrix<T>& a_mat)
  {
    size_t N = a_mat.compute_storage_size();
    out << "(";
    for(int i = 0; i < N; ++i) {
      out << "( ";
      for(int j = 0; j < N; ++j)
        out << a_mat(i, j) << ", ";
      out << "), ";
    }
    out << "), ";
    return out;
  }
};  // namespace libUncertainty

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
inline correlation_store<double>& get_global_correlation_store()
{
  static correlation_store<double> store;
  return store;
}
}  // namespace libUncertainty
