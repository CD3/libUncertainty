#pragma once
#include <algorithm>
#include <iterator>
#include <numeric>

#include "./utils.hpp"

/** @file statistics.hpp
 * @brief
 * @author C.D. Clark III
 * @date 03/15/22
 */

namespace libUncertainty
{
/**
 * @brief Computes the average of a sample.
 *
 * Works with Boost.Units quantities
 *
 * @param begin iterator pointing to beginning of range.
 * @param end iterator pointing to end of range.
 */
template<typename iterator>
auto average(iterator begin, iterator end)
    -> decltype(zero<typename std::iterator_traits<iterator>::value_type>())
{
  using return_type = typename std::iterator_traits<iterator>::value_type;
  auto        N     = std::distance(begin, end);
  return_type sum   = std::accumulate(begin, end, zero<return_type>());

  return sum / N;
}

/**
 * @brief Computes the variance of a sample (the average of the deviations squared).
 *
 * Works with Boost.Units quantities
 *
 * @param begin iterator pointing to beginning of range.
 * @param end iterator pointing to end of range.
 * @param degree_of_freedom_reduce the reduction in the number of degrees of freedom to use. 1 (default) gives the "unbiast estimate" of variance, 0 gives the biast estimate.
 */
template<typename iterator>
auto variance(iterator begin, iterator end, size_t degree_of_freedom_reduce = 1)
    -> decltype(zero<typename std::iterator_traits<iterator>::value_type>() * zero<typename std::iterator_traits<iterator>::value_type>())
{
  using value_type  = typename std::iterator_traits<iterator>::value_type;
  using return_type = decltype(zero<value_type>() * zero<value_type>());
  auto        N     = std::distance(begin, end);
  auto        mu    = average(begin, end);
  return_type sum   = std::accumulate(begin, end, zero<return_type>(), [&mu](return_type sum, value_type x) { return std::move(sum) + (x - mu) * (x - mu); });

  return sum / (N - degree_of_freedom_reduce);
}

/**
 * @brief Computes the standard deviation of a sample (the squre root of the variance).
 *
 * Works with Boost.Units quantities
 *
 * @param begin iterator pointing to beginning of range.
 * @param end iterator pointing to end of range.
 * @param degree_of_freedom_reduce the reduction in the number of degrees of freedom to use. 1 (default) gives the square root of the "unbiast estimate" of variance, 0 gives the biast estimate.
 */
template<typename iterator>
auto standard_deviation(iterator begin, iterator end, size_t degree_of_freedom_reduce = 1)
    -> decltype(sqrt(variance(begin, end, degree_of_freedom_reduce)))
{
  return sqrt(variance(begin, end, degree_of_freedom_reduce));
}

/**
 * @brief Computes the standard error of the mean (the standard deviation divided by the square root of the number of samples).
 *
 * Works with Boost.Units quantities
 *
 * @param begin iterator pointing to beginning of range.
 * @param end iterator pointing to end of range.
 * @param degree_of_freedom_reduce the reduction in the number of degrees of freedom to use. 1 (default) gives the square root of the "unbiast estimate" of variance, 0 gives the biast estimate.
 */
template<typename iterator>
auto standard_error_of_the_mean(iterator begin, iterator end)
    -> decltype(standard_deviation(begin, end, 1))
{
  double N = std::distance(begin, end);
  return standard_deviation(begin, end, 1) / sqrt(N);
}

/**
 * Compute the z-score between two values.
 *
 * z = \frac{|m_1 - m_1|}{\sqrt{\Delta m_1^2 + \Delta m_2^2}}
 */
template<typename T1, typename T2>
double z_score(const T1& a_a, const T2& a_b)

{
  auto a_unc = get_uncertainty(a_a);
  auto b_unc = get_uncertainty(a_b);
  return abs(get_nominal(a_a) - get_nominal(a_b)) / sqrt(a_unc * a_unc + b_unc * b_unc);
}

}  // namespace libUncertainty
