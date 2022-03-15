#pragma once
#include <array>
#include <numeric>

#include "./correlation.hpp"
#include "./uncertain.hpp"
#include "./utils.hpp"

/** @file propagate.hpp
  * @brief 
  * @author C.D. Clark III
  * @date 03/05/22
  */

namespace libUncertainty
{
  /**
   * A class that provides basic error propagation through arbitrary functions.
   */
struct basic_error_propagator {
  template<typename T, size_t N>
  using static_vector = std::array<T, N>;

  /**
   * Propagate error through a function f.
   *
   * DOES NOT HANDLE CORRELATED INPUTS
   */
  template<typename F, typename... Args>
  static auto propagate_error(F a_f, Args... args)
      -> uncertain<decltype(a_f(get_nominal(args)...))>
  {
    // [1] Need to be careful here. It is possible that the difference between
    // two returned values has a different type than a single return value.
    // For example, if the function returns a type representing a quantity
    // with a unit that has an offset (i.e. temperature in celcius: 100 C - 90 C 10 delta_C \ne 10 C)
    static_vector<decltype(a_f(get_nominal(args)...) - a_f(get_upper(args)...)), sizeof...(Args)> deviations;
    auto                                                                                    nominal = _propagate_error(a_f, deviations, std::forward<Args>(args)...);
    auto                                                                                    unc     = sqrt(std::inner_product(deviations.begin() + 1, deviations.end(), deviations.begin() + 1, deviations[0] * deviations[0]));
    uncertain<decltype(a_f(get_nominal(args)...))>                                             ret(nominal, unc);
    return ret;
  }

  /**
   * Propagate error through a function f with a correlations passed in as a matrix.
   */
  template<typename F, typename CorrelationMatrixType, typename... Args>
  static auto propagate_error(F a_f, const CorrelationMatrixType& a_correlation_matrix, Args... args)
      -> uncertain<decltype(a_f(args.nominal()...))>
  {
    // See note [1] above
    static_vector<decltype(a_f(args.nominal()...) - a_f(args.upper()...)), sizeof...(Args)> deviations;
    auto                                                                                    nominal = _propagate_error(a_f, deviations, std::forward<Args>(args)...);
    auto                                                                                    sum     = std::inner_product(deviations.begin() + 1, deviations.end(), deviations.begin() + 1, deviations[0] * deviations[0]);
    for(int i = 0; i < deviations.size(); i++) {
      for(int j = i + 1; j < deviations.size(); j++) {
        sum += 2 * a_correlation_matrix(i, j) * deviations[i] * deviations[j];
      }
    }
    auto                                        unc = sqrt(sum);
    uncertain<decltype(a_f(args.nominal()...))> ret(nominal, unc);
    return ret;
  }

  /**
   * Propagate error through a function f and returns the result with correlations.
   *
   * DOES NOT HANDLE CORRELATED INPUTS
   */
  template<typename F, typename... Args>
  static auto propagate_error_and_correlation(F a_f, Args... args)
      -> add_correlation_coefficient_array<uncertain<decltype(a_f(args.nominal()...))>, double>
  {
    // See note [1] above
    static_vector<decltype(a_f(args.nominal()...) - a_f(args.upper()...)), sizeof...(Args)> deviations;
    auto                                                                                    nominal = _propagate_error(a_f, deviations, std::forward<Args>(args)...);
    auto                                                                                    unc     = sqrt(std::inner_product(deviations.begin() + 1, deviations.end(), deviations.begin() + 1, deviations[0] * deviations[0]));
    add_correlation_coefficient_array<uncertain<decltype(a_f(args.nominal()...))>, double>  ret(nominal, unc);
    ret.set_correlation_coefficient_array_size(sizeof...(Args));
    std::transform(deviations.begin(), deviations.end(), ret.get_correlation_coefficients().begin(),
                   [&unc](auto dev) { return dev / unc; });
    return ret;
  }

  /**
   * Propagate error through a function f with correlations and returns the result with correlations.
   */
  template<typename F, typename CorrelationMatrixType, typename... Args>
  static auto propagate_error_and_correlation(F a_f, const CorrelationMatrixType& a_correlation_matrix, Args... args)
      -> add_correlation_coefficient_array<uncertain<decltype(a_f(args.nominal()...))>, double>
  {
    // See note [1] above
    using deviations_type = decltype(a_f(args.nominal()...) - a_f(args.upper()...));
    static_vector<deviations_type, sizeof...(Args)> deviations;
    auto                                            nominal = _propagate_error(a_f, deviations, std::forward<Args>(args)...);
    auto                                            sum     = std::inner_product(deviations.begin() + 1, deviations.end(), deviations.begin() + 1, deviations[0] * deviations[0]);
    for(int i = 0; i < deviations.size(); i++) {
      for(int j = i + 1; j < deviations.size(); j++) {
        sum += 2 * a_correlation_matrix(i, j) * deviations[i] * deviations[j];
      }
    }
    auto                                                                                   unc = sqrt(sum);
    add_correlation_coefficient_array<uncertain<decltype(a_f(args.nominal()...))>, double> ret(nominal, unc);
    ret.set_correlation_coefficient_array_size(sizeof...(Args));
    for(int i = 0; i < sizeof...(Args); ++i) {
      auto sum = deviations[i];
      for(int j = 0; j < sizeof...(Args); ++j) {
        if(i != j) {
          sum += a_correlation_matrix(i, j) * deviations[j];
        }
      }
      ret.get_correlation_coefficients()[i] = sum / unc;
    }
    return ret;
  }

  /**
   * Propagate error through a function f with correlations using a correlation store.
   */
  template<typename F, typename T, typename... Args>
  static auto propagate_error(F a_f, correlation_store<T>& a_correlation_store, Args... args)
      -> add_id<uncertain<decltype(a_f(args.nominal()...))>>
  {
    // See note [1] above
    using deviations_type = decltype(a_f(args.nominal()...) - a_f(args.upper()...));
    using id_type         = decltype(get_uniq_id());

    static_vector<deviations_type, sizeof...(Args)> deviations;
    static_vector<id_type, sizeof...(Args)>         ids{get_id(args)...};

    auto nominal = _propagate_error(a_f, deviations, std::forward<Args>(args)...);
    
    // compute the uncertainty
    // start with uncorrelated terms
    auto sum     = std::inner_product(deviations.begin() + 1, deviations.end(), deviations.begin() + 1, deviations[0] * deviations[0]);
    // add correlation terms
    for(int i = 0; i < deviations.size(); i++) {
      if(ids[i] == 0) {
        continue;
      }
      for(int j = i + 1; j < deviations.size(); j++) {
        if(ids[j] == 0) {
          continue;
        }
        sum += 2 * a_correlation_store.get_with_ids(ids[i], ids[j]) * deviations[i] * deviations[j];
      }
    }
    auto                                                unc = sqrt(sum);

    // return value
    add_id<uncertain<decltype(a_f(args.nominal()...))>> ret(nominal, unc);

    // compute correlation coefficients for the return value and add them to the store
    for(int i = 0; i < sizeof...(Args); ++i) {
      auto sum = deviations[i];
      for(int j = 0; j < sizeof...(Args); ++j) {
        if(i != j) {
          sum += a_correlation_store.get_with_ids(ids[i],ids[j])*deviations[j];
        }
      }
      a_correlation_store.set_with_ids( ret.get_id(), ids[i], sum/unc);
    }

    return ret;
  }

 private:
// BEGIN GENERATED CODE
// this code was generated using the generate_basic_error_propagator_propagate_error_templates.py script
// to delete this code, you can run (in vim) :g/^\s*\/\/ BEGIN GENERATED CODE/,/^\s*\/\/ END GENERATED CODE/ d
// to reinsert it, you can run (in vim) :.! python scripts/generate_basic_error_propagator_propagate_error_templates.py
// I'm sorry this is so old-school, but your debugger will thank me...
template<typename F, typename T,
typename A0
>
static auto _propagate_error(F a_f,
static_vector<T,1>& a_deviations,
const A0& a_a0
)
{
auto nominal = a_f( get_nominal(a_a0));
if( is_uncertain(a_a0) ) {
a_deviations[0]= a_f( get_upper(a_a0)) - nominal;
}else{
a_deviations[0] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
return nominal;
}

template<typename F, typename T,
typename A0,
typename A1
>
static auto _propagate_error(F a_f,
static_vector<T,2>& a_deviations,
const A0& a_a0,
const A1& a_a1
)
{
auto nominal = a_f( get_nominal(a_a0), get_nominal(a_a1));
if( is_uncertain(a_a0) ) {
a_deviations[0]= a_f( get_upper(a_a0), get_nominal(a_a1)) - nominal;
}else{
a_deviations[0] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a1) ) {
a_deviations[1]= a_f( get_nominal(a_a0), get_upper(a_a1)) - nominal;
}else{
a_deviations[1] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
return nominal;
}

template<typename F, typename T,
typename A0,
typename A1,
typename A2
>
static auto _propagate_error(F a_f,
static_vector<T,3>& a_deviations,
const A0& a_a0,
const A1& a_a1,
const A2& a_a2
)
{
auto nominal = a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2));
if( is_uncertain(a_a0) ) {
a_deviations[0]= a_f( get_upper(a_a0), get_nominal(a_a1), get_nominal(a_a2)) - nominal;
}else{
a_deviations[0] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a1) ) {
a_deviations[1]= a_f( get_nominal(a_a0), get_upper(a_a1), get_nominal(a_a2)) - nominal;
}else{
a_deviations[1] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a2) ) {
a_deviations[2]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_upper(a_a2)) - nominal;
}else{
a_deviations[2] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
return nominal;
}

template<typename F, typename T,
typename A0,
typename A1,
typename A2,
typename A3
>
static auto _propagate_error(F a_f,
static_vector<T,4>& a_deviations,
const A0& a_a0,
const A1& a_a1,
const A2& a_a2,
const A3& a_a3
)
{
auto nominal = a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3));
if( is_uncertain(a_a0) ) {
a_deviations[0]= a_f( get_upper(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3)) - nominal;
}else{
a_deviations[0] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a1) ) {
a_deviations[1]= a_f( get_nominal(a_a0), get_upper(a_a1), get_nominal(a_a2), get_nominal(a_a3)) - nominal;
}else{
a_deviations[1] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a2) ) {
a_deviations[2]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_upper(a_a2), get_nominal(a_a3)) - nominal;
}else{
a_deviations[2] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a3) ) {
a_deviations[3]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_upper(a_a3)) - nominal;
}else{
a_deviations[3] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
return nominal;
}

template<typename F, typename T,
typename A0,
typename A1,
typename A2,
typename A3,
typename A4
>
static auto _propagate_error(F a_f,
static_vector<T,5>& a_deviations,
const A0& a_a0,
const A1& a_a1,
const A2& a_a2,
const A3& a_a3,
const A4& a_a4
)
{
auto nominal = a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4));
if( is_uncertain(a_a0) ) {
a_deviations[0]= a_f( get_upper(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4)) - nominal;
}else{
a_deviations[0] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a1) ) {
a_deviations[1]= a_f( get_nominal(a_a0), get_upper(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4)) - nominal;
}else{
a_deviations[1] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a2) ) {
a_deviations[2]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_upper(a_a2), get_nominal(a_a3), get_nominal(a_a4)) - nominal;
}else{
a_deviations[2] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a3) ) {
a_deviations[3]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_upper(a_a3), get_nominal(a_a4)) - nominal;
}else{
a_deviations[3] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a4) ) {
a_deviations[4]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_upper(a_a4)) - nominal;
}else{
a_deviations[4] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
return nominal;
}

template<typename F, typename T,
typename A0,
typename A1,
typename A2,
typename A3,
typename A4,
typename A5
>
static auto _propagate_error(F a_f,
static_vector<T,6>& a_deviations,
const A0& a_a0,
const A1& a_a1,
const A2& a_a2,
const A3& a_a3,
const A4& a_a4,
const A5& a_a5
)
{
auto nominal = a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5));
if( is_uncertain(a_a0) ) {
a_deviations[0]= a_f( get_upper(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5)) - nominal;
}else{
a_deviations[0] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a1) ) {
a_deviations[1]= a_f( get_nominal(a_a0), get_upper(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5)) - nominal;
}else{
a_deviations[1] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a2) ) {
a_deviations[2]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_upper(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5)) - nominal;
}else{
a_deviations[2] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a3) ) {
a_deviations[3]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_upper(a_a3), get_nominal(a_a4), get_nominal(a_a5)) - nominal;
}else{
a_deviations[3] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a4) ) {
a_deviations[4]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_upper(a_a4), get_nominal(a_a5)) - nominal;
}else{
a_deviations[4] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a5) ) {
a_deviations[5]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_upper(a_a5)) - nominal;
}else{
a_deviations[5] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
return nominal;
}

template<typename F, typename T,
typename A0,
typename A1,
typename A2,
typename A3,
typename A4,
typename A5,
typename A6
>
static auto _propagate_error(F a_f,
static_vector<T,7>& a_deviations,
const A0& a_a0,
const A1& a_a1,
const A2& a_a2,
const A3& a_a3,
const A4& a_a4,
const A5& a_a5,
const A6& a_a6
)
{
auto nominal = a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6));
if( is_uncertain(a_a0) ) {
a_deviations[0]= a_f( get_upper(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6)) - nominal;
}else{
a_deviations[0] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a1) ) {
a_deviations[1]= a_f( get_nominal(a_a0), get_upper(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6)) - nominal;
}else{
a_deviations[1] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a2) ) {
a_deviations[2]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_upper(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6)) - nominal;
}else{
a_deviations[2] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a3) ) {
a_deviations[3]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_upper(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6)) - nominal;
}else{
a_deviations[3] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a4) ) {
a_deviations[4]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_upper(a_a4), get_nominal(a_a5), get_nominal(a_a6)) - nominal;
}else{
a_deviations[4] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a5) ) {
a_deviations[5]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_upper(a_a5), get_nominal(a_a6)) - nominal;
}else{
a_deviations[5] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a6) ) {
a_deviations[6]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_upper(a_a6)) - nominal;
}else{
a_deviations[6] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
return nominal;
}

template<typename F, typename T,
typename A0,
typename A1,
typename A2,
typename A3,
typename A4,
typename A5,
typename A6,
typename A7
>
static auto _propagate_error(F a_f,
static_vector<T,8>& a_deviations,
const A0& a_a0,
const A1& a_a1,
const A2& a_a2,
const A3& a_a3,
const A4& a_a4,
const A5& a_a5,
const A6& a_a6,
const A7& a_a7
)
{
auto nominal = a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7));
if( is_uncertain(a_a0) ) {
a_deviations[0]= a_f( get_upper(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7)) - nominal;
}else{
a_deviations[0] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a1) ) {
a_deviations[1]= a_f( get_nominal(a_a0), get_upper(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7)) - nominal;
}else{
a_deviations[1] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a2) ) {
a_deviations[2]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_upper(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7)) - nominal;
}else{
a_deviations[2] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a3) ) {
a_deviations[3]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_upper(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7)) - nominal;
}else{
a_deviations[3] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a4) ) {
a_deviations[4]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_upper(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7)) - nominal;
}else{
a_deviations[4] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a5) ) {
a_deviations[5]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_upper(a_a5), get_nominal(a_a6), get_nominal(a_a7)) - nominal;
}else{
a_deviations[5] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a6) ) {
a_deviations[6]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_upper(a_a6), get_nominal(a_a7)) - nominal;
}else{
a_deviations[6] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a7) ) {
a_deviations[7]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_upper(a_a7)) - nominal;
}else{
a_deviations[7] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
return nominal;
}

template<typename F, typename T,
typename A0,
typename A1,
typename A2,
typename A3,
typename A4,
typename A5,
typename A6,
typename A7,
typename A8
>
static auto _propagate_error(F a_f,
static_vector<T,9>& a_deviations,
const A0& a_a0,
const A1& a_a1,
const A2& a_a2,
const A3& a_a3,
const A4& a_a4,
const A5& a_a5,
const A6& a_a6,
const A7& a_a7,
const A8& a_a8
)
{
auto nominal = a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8));
if( is_uncertain(a_a0) ) {
a_deviations[0]= a_f( get_upper(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8)) - nominal;
}else{
a_deviations[0] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a1) ) {
a_deviations[1]= a_f( get_nominal(a_a0), get_upper(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8)) - nominal;
}else{
a_deviations[1] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a2) ) {
a_deviations[2]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_upper(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8)) - nominal;
}else{
a_deviations[2] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a3) ) {
a_deviations[3]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_upper(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8)) - nominal;
}else{
a_deviations[3] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a4) ) {
a_deviations[4]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_upper(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8)) - nominal;
}else{
a_deviations[4] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a5) ) {
a_deviations[5]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_upper(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8)) - nominal;
}else{
a_deviations[5] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a6) ) {
a_deviations[6]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_upper(a_a6), get_nominal(a_a7), get_nominal(a_a8)) - nominal;
}else{
a_deviations[6] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a7) ) {
a_deviations[7]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_upper(a_a7), get_nominal(a_a8)) - nominal;
}else{
a_deviations[7] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a8) ) {
a_deviations[8]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_upper(a_a8)) - nominal;
}else{
a_deviations[8] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
return nominal;
}

template<typename F, typename T,
typename A0,
typename A1,
typename A2,
typename A3,
typename A4,
typename A5,
typename A6,
typename A7,
typename A8,
typename A9
>
static auto _propagate_error(F a_f,
static_vector<T,10>& a_deviations,
const A0& a_a0,
const A1& a_a1,
const A2& a_a2,
const A3& a_a3,
const A4& a_a4,
const A5& a_a5,
const A6& a_a6,
const A7& a_a7,
const A8& a_a8,
const A9& a_a9
)
{
auto nominal = a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9));
if( is_uncertain(a_a0) ) {
a_deviations[0]= a_f( get_upper(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9)) - nominal;
}else{
a_deviations[0] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a1) ) {
a_deviations[1]= a_f( get_nominal(a_a0), get_upper(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9)) - nominal;
}else{
a_deviations[1] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a2) ) {
a_deviations[2]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_upper(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9)) - nominal;
}else{
a_deviations[2] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a3) ) {
a_deviations[3]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_upper(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9)) - nominal;
}else{
a_deviations[3] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a4) ) {
a_deviations[4]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_upper(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9)) - nominal;
}else{
a_deviations[4] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a5) ) {
a_deviations[5]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_upper(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9)) - nominal;
}else{
a_deviations[5] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a6) ) {
a_deviations[6]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_upper(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9)) - nominal;
}else{
a_deviations[6] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a7) ) {
a_deviations[7]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_upper(a_a7), get_nominal(a_a8), get_nominal(a_a9)) - nominal;
}else{
a_deviations[7] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a8) ) {
a_deviations[8]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_upper(a_a8), get_nominal(a_a9)) - nominal;
}else{
a_deviations[8] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a9) ) {
a_deviations[9]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_upper(a_a9)) - nominal;
}else{
a_deviations[9] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
return nominal;
}

template<typename F, typename T,
typename A0,
typename A1,
typename A2,
typename A3,
typename A4,
typename A5,
typename A6,
typename A7,
typename A8,
typename A9,
typename A10
>
static auto _propagate_error(F a_f,
static_vector<T,11>& a_deviations,
const A0& a_a0,
const A1& a_a1,
const A2& a_a2,
const A3& a_a3,
const A4& a_a4,
const A5& a_a5,
const A6& a_a6,
const A7& a_a7,
const A8& a_a8,
const A9& a_a9,
const A10& a_a10
)
{
auto nominal = a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10));
if( is_uncertain(a_a0) ) {
a_deviations[0]= a_f( get_upper(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10)) - nominal;
}else{
a_deviations[0] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a1) ) {
a_deviations[1]= a_f( get_nominal(a_a0), get_upper(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10)) - nominal;
}else{
a_deviations[1] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a2) ) {
a_deviations[2]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_upper(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10)) - nominal;
}else{
a_deviations[2] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a3) ) {
a_deviations[3]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_upper(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10)) - nominal;
}else{
a_deviations[3] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a4) ) {
a_deviations[4]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_upper(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10)) - nominal;
}else{
a_deviations[4] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a5) ) {
a_deviations[5]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_upper(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10)) - nominal;
}else{
a_deviations[5] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a6) ) {
a_deviations[6]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_upper(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10)) - nominal;
}else{
a_deviations[6] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a7) ) {
a_deviations[7]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_upper(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10)) - nominal;
}else{
a_deviations[7] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a8) ) {
a_deviations[8]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_upper(a_a8), get_nominal(a_a9), get_nominal(a_a10)) - nominal;
}else{
a_deviations[8] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a9) ) {
a_deviations[9]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_upper(a_a9), get_nominal(a_a10)) - nominal;
}else{
a_deviations[9] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a10) ) {
a_deviations[10]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_upper(a_a10)) - nominal;
}else{
a_deviations[10] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
return nominal;
}

template<typename F, typename T,
typename A0,
typename A1,
typename A2,
typename A3,
typename A4,
typename A5,
typename A6,
typename A7,
typename A8,
typename A9,
typename A10,
typename A11
>
static auto _propagate_error(F a_f,
static_vector<T,12>& a_deviations,
const A0& a_a0,
const A1& a_a1,
const A2& a_a2,
const A3& a_a3,
const A4& a_a4,
const A5& a_a5,
const A6& a_a6,
const A7& a_a7,
const A8& a_a8,
const A9& a_a9,
const A10& a_a10,
const A11& a_a11
)
{
auto nominal = a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11));
if( is_uncertain(a_a0) ) {
a_deviations[0]= a_f( get_upper(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11)) - nominal;
}else{
a_deviations[0] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a1) ) {
a_deviations[1]= a_f( get_nominal(a_a0), get_upper(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11)) - nominal;
}else{
a_deviations[1] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a2) ) {
a_deviations[2]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_upper(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11)) - nominal;
}else{
a_deviations[2] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a3) ) {
a_deviations[3]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_upper(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11)) - nominal;
}else{
a_deviations[3] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a4) ) {
a_deviations[4]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_upper(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11)) - nominal;
}else{
a_deviations[4] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a5) ) {
a_deviations[5]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_upper(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11)) - nominal;
}else{
a_deviations[5] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a6) ) {
a_deviations[6]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_upper(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11)) - nominal;
}else{
a_deviations[6] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a7) ) {
a_deviations[7]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_upper(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11)) - nominal;
}else{
a_deviations[7] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a8) ) {
a_deviations[8]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_upper(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11)) - nominal;
}else{
a_deviations[8] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a9) ) {
a_deviations[9]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_upper(a_a9), get_nominal(a_a10), get_nominal(a_a11)) - nominal;
}else{
a_deviations[9] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a10) ) {
a_deviations[10]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_upper(a_a10), get_nominal(a_a11)) - nominal;
}else{
a_deviations[10] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a11) ) {
a_deviations[11]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_upper(a_a11)) - nominal;
}else{
a_deviations[11] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
return nominal;
}

template<typename F, typename T,
typename A0,
typename A1,
typename A2,
typename A3,
typename A4,
typename A5,
typename A6,
typename A7,
typename A8,
typename A9,
typename A10,
typename A11,
typename A12
>
static auto _propagate_error(F a_f,
static_vector<T,13>& a_deviations,
const A0& a_a0,
const A1& a_a1,
const A2& a_a2,
const A3& a_a3,
const A4& a_a4,
const A5& a_a5,
const A6& a_a6,
const A7& a_a7,
const A8& a_a8,
const A9& a_a9,
const A10& a_a10,
const A11& a_a11,
const A12& a_a12
)
{
auto nominal = a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12));
if( is_uncertain(a_a0) ) {
a_deviations[0]= a_f( get_upper(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12)) - nominal;
}else{
a_deviations[0] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a1) ) {
a_deviations[1]= a_f( get_nominal(a_a0), get_upper(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12)) - nominal;
}else{
a_deviations[1] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a2) ) {
a_deviations[2]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_upper(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12)) - nominal;
}else{
a_deviations[2] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a3) ) {
a_deviations[3]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_upper(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12)) - nominal;
}else{
a_deviations[3] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a4) ) {
a_deviations[4]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_upper(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12)) - nominal;
}else{
a_deviations[4] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a5) ) {
a_deviations[5]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_upper(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12)) - nominal;
}else{
a_deviations[5] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a6) ) {
a_deviations[6]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_upper(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12)) - nominal;
}else{
a_deviations[6] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a7) ) {
a_deviations[7]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_upper(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12)) - nominal;
}else{
a_deviations[7] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a8) ) {
a_deviations[8]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_upper(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12)) - nominal;
}else{
a_deviations[8] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a9) ) {
a_deviations[9]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_upper(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12)) - nominal;
}else{
a_deviations[9] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a10) ) {
a_deviations[10]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_upper(a_a10), get_nominal(a_a11), get_nominal(a_a12)) - nominal;
}else{
a_deviations[10] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a11) ) {
a_deviations[11]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_upper(a_a11), get_nominal(a_a12)) - nominal;
}else{
a_deviations[11] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a12) ) {
a_deviations[12]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_upper(a_a12)) - nominal;
}else{
a_deviations[12] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
return nominal;
}

template<typename F, typename T,
typename A0,
typename A1,
typename A2,
typename A3,
typename A4,
typename A5,
typename A6,
typename A7,
typename A8,
typename A9,
typename A10,
typename A11,
typename A12,
typename A13
>
static auto _propagate_error(F a_f,
static_vector<T,14>& a_deviations,
const A0& a_a0,
const A1& a_a1,
const A2& a_a2,
const A3& a_a3,
const A4& a_a4,
const A5& a_a5,
const A6& a_a6,
const A7& a_a7,
const A8& a_a8,
const A9& a_a9,
const A10& a_a10,
const A11& a_a11,
const A12& a_a12,
const A13& a_a13
)
{
auto nominal = a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13));
if( is_uncertain(a_a0) ) {
a_deviations[0]= a_f( get_upper(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13)) - nominal;
}else{
a_deviations[0] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a1) ) {
a_deviations[1]= a_f( get_nominal(a_a0), get_upper(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13)) - nominal;
}else{
a_deviations[1] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a2) ) {
a_deviations[2]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_upper(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13)) - nominal;
}else{
a_deviations[2] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a3) ) {
a_deviations[3]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_upper(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13)) - nominal;
}else{
a_deviations[3] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a4) ) {
a_deviations[4]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_upper(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13)) - nominal;
}else{
a_deviations[4] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a5) ) {
a_deviations[5]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_upper(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13)) - nominal;
}else{
a_deviations[5] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a6) ) {
a_deviations[6]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_upper(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13)) - nominal;
}else{
a_deviations[6] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a7) ) {
a_deviations[7]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_upper(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13)) - nominal;
}else{
a_deviations[7] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a8) ) {
a_deviations[8]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_upper(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13)) - nominal;
}else{
a_deviations[8] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a9) ) {
a_deviations[9]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_upper(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13)) - nominal;
}else{
a_deviations[9] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a10) ) {
a_deviations[10]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_upper(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13)) - nominal;
}else{
a_deviations[10] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a11) ) {
a_deviations[11]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_upper(a_a11), get_nominal(a_a12), get_nominal(a_a13)) - nominal;
}else{
a_deviations[11] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a12) ) {
a_deviations[12]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_upper(a_a12), get_nominal(a_a13)) - nominal;
}else{
a_deviations[12] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a13) ) {
a_deviations[13]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_upper(a_a13)) - nominal;
}else{
a_deviations[13] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
return nominal;
}

template<typename F, typename T,
typename A0,
typename A1,
typename A2,
typename A3,
typename A4,
typename A5,
typename A6,
typename A7,
typename A8,
typename A9,
typename A10,
typename A11,
typename A12,
typename A13,
typename A14
>
static auto _propagate_error(F a_f,
static_vector<T,15>& a_deviations,
const A0& a_a0,
const A1& a_a1,
const A2& a_a2,
const A3& a_a3,
const A4& a_a4,
const A5& a_a5,
const A6& a_a6,
const A7& a_a7,
const A8& a_a8,
const A9& a_a9,
const A10& a_a10,
const A11& a_a11,
const A12& a_a12,
const A13& a_a13,
const A14& a_a14
)
{
auto nominal = a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14));
if( is_uncertain(a_a0) ) {
a_deviations[0]= a_f( get_upper(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14)) - nominal;
}else{
a_deviations[0] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a1) ) {
a_deviations[1]= a_f( get_nominal(a_a0), get_upper(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14)) - nominal;
}else{
a_deviations[1] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a2) ) {
a_deviations[2]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_upper(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14)) - nominal;
}else{
a_deviations[2] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a3) ) {
a_deviations[3]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_upper(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14)) - nominal;
}else{
a_deviations[3] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a4) ) {
a_deviations[4]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_upper(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14)) - nominal;
}else{
a_deviations[4] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a5) ) {
a_deviations[5]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_upper(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14)) - nominal;
}else{
a_deviations[5] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a6) ) {
a_deviations[6]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_upper(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14)) - nominal;
}else{
a_deviations[6] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a7) ) {
a_deviations[7]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_upper(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14)) - nominal;
}else{
a_deviations[7] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a8) ) {
a_deviations[8]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_upper(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14)) - nominal;
}else{
a_deviations[8] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a9) ) {
a_deviations[9]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_upper(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14)) - nominal;
}else{
a_deviations[9] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a10) ) {
a_deviations[10]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_upper(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14)) - nominal;
}else{
a_deviations[10] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a11) ) {
a_deviations[11]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_upper(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14)) - nominal;
}else{
a_deviations[11] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a12) ) {
a_deviations[12]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_upper(a_a12), get_nominal(a_a13), get_nominal(a_a14)) - nominal;
}else{
a_deviations[12] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a13) ) {
a_deviations[13]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_upper(a_a13), get_nominal(a_a14)) - nominal;
}else{
a_deviations[13] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a14) ) {
a_deviations[14]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_upper(a_a14)) - nominal;
}else{
a_deviations[14] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
return nominal;
}

template<typename F, typename T,
typename A0,
typename A1,
typename A2,
typename A3,
typename A4,
typename A5,
typename A6,
typename A7,
typename A8,
typename A9,
typename A10,
typename A11,
typename A12,
typename A13,
typename A14,
typename A15
>
static auto _propagate_error(F a_f,
static_vector<T,16>& a_deviations,
const A0& a_a0,
const A1& a_a1,
const A2& a_a2,
const A3& a_a3,
const A4& a_a4,
const A5& a_a5,
const A6& a_a6,
const A7& a_a7,
const A8& a_a8,
const A9& a_a9,
const A10& a_a10,
const A11& a_a11,
const A12& a_a12,
const A13& a_a13,
const A14& a_a14,
const A15& a_a15
)
{
auto nominal = a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15));
if( is_uncertain(a_a0) ) {
a_deviations[0]= a_f( get_upper(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15)) - nominal;
}else{
a_deviations[0] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a1) ) {
a_deviations[1]= a_f( get_nominal(a_a0), get_upper(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15)) - nominal;
}else{
a_deviations[1] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a2) ) {
a_deviations[2]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_upper(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15)) - nominal;
}else{
a_deviations[2] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a3) ) {
a_deviations[3]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_upper(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15)) - nominal;
}else{
a_deviations[3] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a4) ) {
a_deviations[4]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_upper(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15)) - nominal;
}else{
a_deviations[4] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a5) ) {
a_deviations[5]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_upper(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15)) - nominal;
}else{
a_deviations[5] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a6) ) {
a_deviations[6]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_upper(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15)) - nominal;
}else{
a_deviations[6] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a7) ) {
a_deviations[7]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_upper(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15)) - nominal;
}else{
a_deviations[7] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a8) ) {
a_deviations[8]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_upper(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15)) - nominal;
}else{
a_deviations[8] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a9) ) {
a_deviations[9]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_upper(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15)) - nominal;
}else{
a_deviations[9] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a10) ) {
a_deviations[10]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_upper(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15)) - nominal;
}else{
a_deviations[10] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a11) ) {
a_deviations[11]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_upper(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15)) - nominal;
}else{
a_deviations[11] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a12) ) {
a_deviations[12]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_upper(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15)) - nominal;
}else{
a_deviations[12] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a13) ) {
a_deviations[13]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_upper(a_a13), get_nominal(a_a14), get_nominal(a_a15)) - nominal;
}else{
a_deviations[13] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a14) ) {
a_deviations[14]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_upper(a_a14), get_nominal(a_a15)) - nominal;
}else{
a_deviations[14] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a15) ) {
a_deviations[15]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_upper(a_a15)) - nominal;
}else{
a_deviations[15] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
return nominal;
}

template<typename F, typename T,
typename A0,
typename A1,
typename A2,
typename A3,
typename A4,
typename A5,
typename A6,
typename A7,
typename A8,
typename A9,
typename A10,
typename A11,
typename A12,
typename A13,
typename A14,
typename A15,
typename A16
>
static auto _propagate_error(F a_f,
static_vector<T,17>& a_deviations,
const A0& a_a0,
const A1& a_a1,
const A2& a_a2,
const A3& a_a3,
const A4& a_a4,
const A5& a_a5,
const A6& a_a6,
const A7& a_a7,
const A8& a_a8,
const A9& a_a9,
const A10& a_a10,
const A11& a_a11,
const A12& a_a12,
const A13& a_a13,
const A14& a_a14,
const A15& a_a15,
const A16& a_a16
)
{
auto nominal = a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_nominal(a_a16));
if( is_uncertain(a_a0) ) {
a_deviations[0]= a_f( get_upper(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_nominal(a_a16)) - nominal;
}else{
a_deviations[0] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a1) ) {
a_deviations[1]= a_f( get_nominal(a_a0), get_upper(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_nominal(a_a16)) - nominal;
}else{
a_deviations[1] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a2) ) {
a_deviations[2]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_upper(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_nominal(a_a16)) - nominal;
}else{
a_deviations[2] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a3) ) {
a_deviations[3]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_upper(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_nominal(a_a16)) - nominal;
}else{
a_deviations[3] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a4) ) {
a_deviations[4]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_upper(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_nominal(a_a16)) - nominal;
}else{
a_deviations[4] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a5) ) {
a_deviations[5]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_upper(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_nominal(a_a16)) - nominal;
}else{
a_deviations[5] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a6) ) {
a_deviations[6]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_upper(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_nominal(a_a16)) - nominal;
}else{
a_deviations[6] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a7) ) {
a_deviations[7]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_upper(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_nominal(a_a16)) - nominal;
}else{
a_deviations[7] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a8) ) {
a_deviations[8]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_upper(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_nominal(a_a16)) - nominal;
}else{
a_deviations[8] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a9) ) {
a_deviations[9]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_upper(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_nominal(a_a16)) - nominal;
}else{
a_deviations[9] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a10) ) {
a_deviations[10]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_upper(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_nominal(a_a16)) - nominal;
}else{
a_deviations[10] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a11) ) {
a_deviations[11]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_upper(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_nominal(a_a16)) - nominal;
}else{
a_deviations[11] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a12) ) {
a_deviations[12]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_upper(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_nominal(a_a16)) - nominal;
}else{
a_deviations[12] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a13) ) {
a_deviations[13]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_upper(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_nominal(a_a16)) - nominal;
}else{
a_deviations[13] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a14) ) {
a_deviations[14]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_upper(a_a14), get_nominal(a_a15), get_nominal(a_a16)) - nominal;
}else{
a_deviations[14] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a15) ) {
a_deviations[15]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_upper(a_a15), get_nominal(a_a16)) - nominal;
}else{
a_deviations[15] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a16) ) {
a_deviations[16]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_upper(a_a16)) - nominal;
}else{
a_deviations[16] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
return nominal;
}

template<typename F, typename T,
typename A0,
typename A1,
typename A2,
typename A3,
typename A4,
typename A5,
typename A6,
typename A7,
typename A8,
typename A9,
typename A10,
typename A11,
typename A12,
typename A13,
typename A14,
typename A15,
typename A16,
typename A17
>
static auto _propagate_error(F a_f,
static_vector<T,18>& a_deviations,
const A0& a_a0,
const A1& a_a1,
const A2& a_a2,
const A3& a_a3,
const A4& a_a4,
const A5& a_a5,
const A6& a_a6,
const A7& a_a7,
const A8& a_a8,
const A9& a_a9,
const A10& a_a10,
const A11& a_a11,
const A12& a_a12,
const A13& a_a13,
const A14& a_a14,
const A15& a_a15,
const A16& a_a16,
const A17& a_a17
)
{
auto nominal = a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_nominal(a_a16), get_nominal(a_a17));
if( is_uncertain(a_a0) ) {
a_deviations[0]= a_f( get_upper(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_nominal(a_a16), get_nominal(a_a17)) - nominal;
}else{
a_deviations[0] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a1) ) {
a_deviations[1]= a_f( get_nominal(a_a0), get_upper(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_nominal(a_a16), get_nominal(a_a17)) - nominal;
}else{
a_deviations[1] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a2) ) {
a_deviations[2]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_upper(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_nominal(a_a16), get_nominal(a_a17)) - nominal;
}else{
a_deviations[2] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a3) ) {
a_deviations[3]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_upper(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_nominal(a_a16), get_nominal(a_a17)) - nominal;
}else{
a_deviations[3] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a4) ) {
a_deviations[4]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_upper(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_nominal(a_a16), get_nominal(a_a17)) - nominal;
}else{
a_deviations[4] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a5) ) {
a_deviations[5]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_upper(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_nominal(a_a16), get_nominal(a_a17)) - nominal;
}else{
a_deviations[5] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a6) ) {
a_deviations[6]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_upper(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_nominal(a_a16), get_nominal(a_a17)) - nominal;
}else{
a_deviations[6] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a7) ) {
a_deviations[7]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_upper(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_nominal(a_a16), get_nominal(a_a17)) - nominal;
}else{
a_deviations[7] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a8) ) {
a_deviations[8]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_upper(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_nominal(a_a16), get_nominal(a_a17)) - nominal;
}else{
a_deviations[8] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a9) ) {
a_deviations[9]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_upper(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_nominal(a_a16), get_nominal(a_a17)) - nominal;
}else{
a_deviations[9] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a10) ) {
a_deviations[10]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_upper(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_nominal(a_a16), get_nominal(a_a17)) - nominal;
}else{
a_deviations[10] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a11) ) {
a_deviations[11]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_upper(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_nominal(a_a16), get_nominal(a_a17)) - nominal;
}else{
a_deviations[11] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a12) ) {
a_deviations[12]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_upper(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_nominal(a_a16), get_nominal(a_a17)) - nominal;
}else{
a_deviations[12] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a13) ) {
a_deviations[13]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_upper(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_nominal(a_a16), get_nominal(a_a17)) - nominal;
}else{
a_deviations[13] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a14) ) {
a_deviations[14]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_upper(a_a14), get_nominal(a_a15), get_nominal(a_a16), get_nominal(a_a17)) - nominal;
}else{
a_deviations[14] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a15) ) {
a_deviations[15]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_upper(a_a15), get_nominal(a_a16), get_nominal(a_a17)) - nominal;
}else{
a_deviations[15] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a16) ) {
a_deviations[16]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_upper(a_a16), get_nominal(a_a17)) - nominal;
}else{
a_deviations[16] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a17) ) {
a_deviations[17]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_nominal(a_a16), get_upper(a_a17)) - nominal;
}else{
a_deviations[17] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
return nominal;
}

template<typename F, typename T,
typename A0,
typename A1,
typename A2,
typename A3,
typename A4,
typename A5,
typename A6,
typename A7,
typename A8,
typename A9,
typename A10,
typename A11,
typename A12,
typename A13,
typename A14,
typename A15,
typename A16,
typename A17,
typename A18
>
static auto _propagate_error(F a_f,
static_vector<T,19>& a_deviations,
const A0& a_a0,
const A1& a_a1,
const A2& a_a2,
const A3& a_a3,
const A4& a_a4,
const A5& a_a5,
const A6& a_a6,
const A7& a_a7,
const A8& a_a8,
const A9& a_a9,
const A10& a_a10,
const A11& a_a11,
const A12& a_a12,
const A13& a_a13,
const A14& a_a14,
const A15& a_a15,
const A16& a_a16,
const A17& a_a17,
const A18& a_a18
)
{
auto nominal = a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_nominal(a_a16), get_nominal(a_a17), get_nominal(a_a18));
if( is_uncertain(a_a0) ) {
a_deviations[0]= a_f( get_upper(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_nominal(a_a16), get_nominal(a_a17), get_nominal(a_a18)) - nominal;
}else{
a_deviations[0] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a1) ) {
a_deviations[1]= a_f( get_nominal(a_a0), get_upper(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_nominal(a_a16), get_nominal(a_a17), get_nominal(a_a18)) - nominal;
}else{
a_deviations[1] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a2) ) {
a_deviations[2]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_upper(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_nominal(a_a16), get_nominal(a_a17), get_nominal(a_a18)) - nominal;
}else{
a_deviations[2] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a3) ) {
a_deviations[3]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_upper(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_nominal(a_a16), get_nominal(a_a17), get_nominal(a_a18)) - nominal;
}else{
a_deviations[3] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a4) ) {
a_deviations[4]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_upper(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_nominal(a_a16), get_nominal(a_a17), get_nominal(a_a18)) - nominal;
}else{
a_deviations[4] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a5) ) {
a_deviations[5]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_upper(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_nominal(a_a16), get_nominal(a_a17), get_nominal(a_a18)) - nominal;
}else{
a_deviations[5] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a6) ) {
a_deviations[6]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_upper(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_nominal(a_a16), get_nominal(a_a17), get_nominal(a_a18)) - nominal;
}else{
a_deviations[6] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a7) ) {
a_deviations[7]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_upper(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_nominal(a_a16), get_nominal(a_a17), get_nominal(a_a18)) - nominal;
}else{
a_deviations[7] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a8) ) {
a_deviations[8]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_upper(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_nominal(a_a16), get_nominal(a_a17), get_nominal(a_a18)) - nominal;
}else{
a_deviations[8] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a9) ) {
a_deviations[9]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_upper(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_nominal(a_a16), get_nominal(a_a17), get_nominal(a_a18)) - nominal;
}else{
a_deviations[9] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a10) ) {
a_deviations[10]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_upper(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_nominal(a_a16), get_nominal(a_a17), get_nominal(a_a18)) - nominal;
}else{
a_deviations[10] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a11) ) {
a_deviations[11]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_upper(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_nominal(a_a16), get_nominal(a_a17), get_nominal(a_a18)) - nominal;
}else{
a_deviations[11] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a12) ) {
a_deviations[12]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_upper(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_nominal(a_a16), get_nominal(a_a17), get_nominal(a_a18)) - nominal;
}else{
a_deviations[12] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a13) ) {
a_deviations[13]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_upper(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_nominal(a_a16), get_nominal(a_a17), get_nominal(a_a18)) - nominal;
}else{
a_deviations[13] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a14) ) {
a_deviations[14]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_upper(a_a14), get_nominal(a_a15), get_nominal(a_a16), get_nominal(a_a17), get_nominal(a_a18)) - nominal;
}else{
a_deviations[14] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a15) ) {
a_deviations[15]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_upper(a_a15), get_nominal(a_a16), get_nominal(a_a17), get_nominal(a_a18)) - nominal;
}else{
a_deviations[15] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a16) ) {
a_deviations[16]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_upper(a_a16), get_nominal(a_a17), get_nominal(a_a18)) - nominal;
}else{
a_deviations[16] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a17) ) {
a_deviations[17]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_nominal(a_a16), get_upper(a_a17), get_nominal(a_a18)) - nominal;
}else{
a_deviations[17] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a18) ) {
a_deviations[18]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_nominal(a_a16), get_nominal(a_a17), get_upper(a_a18)) - nominal;
}else{
a_deviations[18] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
return nominal;
}

template<typename F, typename T,
typename A0,
typename A1,
typename A2,
typename A3,
typename A4,
typename A5,
typename A6,
typename A7,
typename A8,
typename A9,
typename A10,
typename A11,
typename A12,
typename A13,
typename A14,
typename A15,
typename A16,
typename A17,
typename A18,
typename A19
>
static auto _propagate_error(F a_f,
static_vector<T,20>& a_deviations,
const A0& a_a0,
const A1& a_a1,
const A2& a_a2,
const A3& a_a3,
const A4& a_a4,
const A5& a_a5,
const A6& a_a6,
const A7& a_a7,
const A8& a_a8,
const A9& a_a9,
const A10& a_a10,
const A11& a_a11,
const A12& a_a12,
const A13& a_a13,
const A14& a_a14,
const A15& a_a15,
const A16& a_a16,
const A17& a_a17,
const A18& a_a18,
const A19& a_a19
)
{
auto nominal = a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_nominal(a_a16), get_nominal(a_a17), get_nominal(a_a18), get_nominal(a_a19));
if( is_uncertain(a_a0) ) {
a_deviations[0]= a_f( get_upper(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_nominal(a_a16), get_nominal(a_a17), get_nominal(a_a18), get_nominal(a_a19)) - nominal;
}else{
a_deviations[0] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a1) ) {
a_deviations[1]= a_f( get_nominal(a_a0), get_upper(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_nominal(a_a16), get_nominal(a_a17), get_nominal(a_a18), get_nominal(a_a19)) - nominal;
}else{
a_deviations[1] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a2) ) {
a_deviations[2]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_upper(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_nominal(a_a16), get_nominal(a_a17), get_nominal(a_a18), get_nominal(a_a19)) - nominal;
}else{
a_deviations[2] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a3) ) {
a_deviations[3]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_upper(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_nominal(a_a16), get_nominal(a_a17), get_nominal(a_a18), get_nominal(a_a19)) - nominal;
}else{
a_deviations[3] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a4) ) {
a_deviations[4]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_upper(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_nominal(a_a16), get_nominal(a_a17), get_nominal(a_a18), get_nominal(a_a19)) - nominal;
}else{
a_deviations[4] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a5) ) {
a_deviations[5]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_upper(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_nominal(a_a16), get_nominal(a_a17), get_nominal(a_a18), get_nominal(a_a19)) - nominal;
}else{
a_deviations[5] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a6) ) {
a_deviations[6]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_upper(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_nominal(a_a16), get_nominal(a_a17), get_nominal(a_a18), get_nominal(a_a19)) - nominal;
}else{
a_deviations[6] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a7) ) {
a_deviations[7]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_upper(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_nominal(a_a16), get_nominal(a_a17), get_nominal(a_a18), get_nominal(a_a19)) - nominal;
}else{
a_deviations[7] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a8) ) {
a_deviations[8]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_upper(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_nominal(a_a16), get_nominal(a_a17), get_nominal(a_a18), get_nominal(a_a19)) - nominal;
}else{
a_deviations[8] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a9) ) {
a_deviations[9]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_upper(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_nominal(a_a16), get_nominal(a_a17), get_nominal(a_a18), get_nominal(a_a19)) - nominal;
}else{
a_deviations[9] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a10) ) {
a_deviations[10]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_upper(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_nominal(a_a16), get_nominal(a_a17), get_nominal(a_a18), get_nominal(a_a19)) - nominal;
}else{
a_deviations[10] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a11) ) {
a_deviations[11]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_upper(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_nominal(a_a16), get_nominal(a_a17), get_nominal(a_a18), get_nominal(a_a19)) - nominal;
}else{
a_deviations[11] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a12) ) {
a_deviations[12]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_upper(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_nominal(a_a16), get_nominal(a_a17), get_nominal(a_a18), get_nominal(a_a19)) - nominal;
}else{
a_deviations[12] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a13) ) {
a_deviations[13]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_upper(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_nominal(a_a16), get_nominal(a_a17), get_nominal(a_a18), get_nominal(a_a19)) - nominal;
}else{
a_deviations[13] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a14) ) {
a_deviations[14]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_upper(a_a14), get_nominal(a_a15), get_nominal(a_a16), get_nominal(a_a17), get_nominal(a_a18), get_nominal(a_a19)) - nominal;
}else{
a_deviations[14] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a15) ) {
a_deviations[15]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_upper(a_a15), get_nominal(a_a16), get_nominal(a_a17), get_nominal(a_a18), get_nominal(a_a19)) - nominal;
}else{
a_deviations[15] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a16) ) {
a_deviations[16]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_upper(a_a16), get_nominal(a_a17), get_nominal(a_a18), get_nominal(a_a19)) - nominal;
}else{
a_deviations[16] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a17) ) {
a_deviations[17]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_nominal(a_a16), get_upper(a_a17), get_nominal(a_a18), get_nominal(a_a19)) - nominal;
}else{
a_deviations[17] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a18) ) {
a_deviations[18]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_nominal(a_a16), get_nominal(a_a17), get_upper(a_a18), get_nominal(a_a19)) - nominal;
}else{
a_deviations[18] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
if( is_uncertain(a_a19) ) {
a_deviations[19]= a_f( get_nominal(a_a0), get_nominal(a_a1), get_nominal(a_a2), get_nominal(a_a3), get_nominal(a_a4), get_nominal(a_a5), get_nominal(a_a6), get_nominal(a_a7), get_nominal(a_a8), get_nominal(a_a9), get_nominal(a_a10), get_nominal(a_a11), get_nominal(a_a12), get_nominal(a_a13), get_nominal(a_a14), get_nominal(a_a15), get_nominal(a_a16), get_nominal(a_a17), get_nominal(a_a18), get_upper(a_a19)) - nominal;
}else{
a_deviations[19] = nominal - nominal; // can't just use 0 here, we might be dealing with quantities that have units
}
return nominal;
}

// END GENERATED CODE
};

}  // namespace libUncertainty
