#pragma once
#include <array>
#include <numeric>

#include <boost/numeric/ublas/vector.hpp>

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
  using static_vector = boost::numeric::ublas::fixed_vector<T, N>;

  /**
   * Propagate error through a function f.
   *
   * DOES NOT HANDLE CORRELATED INPUTS
   */
  template<typename F, typename... Args>
  static auto propagate_error(F a_f, Args... args)
      -> uncertain<decltype(a_f(args.nominal()...))>
  {
    // [1] Need to be careful here. It is possible that the difference between
    // two returned values has a different type than a single return value.
    // For example, if the function returns a type representing a quantity
    // with a unit that has an offset (i.e. temperature in celcius: 100 C - 90 C 10 delta_C \ne 10 C)
    static_vector<decltype(a_f(args.nominal()...) - a_f(args.upper()...)), sizeof...(Args)> deviations;
    auto                                                                                    nominal = _propagate_error(a_f, deviations, std::forward<Args>(args)...);
    auto                                                                                    unc     = sqrt(std::inner_product(deviations.begin() + 1, deviations.end(), deviations.begin() + 1, deviations[0] * deviations[0]));
    uncertain<decltype(a_f(args.nominal()...))>                                             ret(nominal, unc);
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
           typename N0, typename U0>
  static auto _propagate_error(F                        a_f,
                               static_vector<T, 1>&     a_deviations,
                               const uncertain<N0, U0>& a_a0)
  {
    auto nominal    = a_f(a_a0.nominal());
    a_deviations[0] = a_f(a_a0.upper()) - nominal;
    return nominal;
  }

  template<typename F, typename T,
           typename N0, typename U0,
           typename N1, typename U1>
  static auto _propagate_error(F                        a_f,
                               static_vector<T, 2>&     a_deviations,
                               const uncertain<N0, U0>& a_a0,
                               const uncertain<N1, U1>& a_a1)
  {
    auto nominal    = a_f(a_a0.nominal(), a_a1.nominal());
    a_deviations[0] = a_f(a_a0.upper(), a_a1.nominal()) - nominal;
    a_deviations[1] = a_f(a_a0.nominal(), a_a1.upper()) - nominal;
    return nominal;
  }

  template<typename F, typename T,
           typename N0, typename U0,
           typename N1, typename U1,
           typename N2, typename U2>
  static auto _propagate_error(F                        a_f,
                               static_vector<T, 3>&     a_deviations,
                               const uncertain<N0, U0>& a_a0,
                               const uncertain<N1, U1>& a_a1,
                               const uncertain<N2, U2>& a_a2)
  {
    auto nominal    = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal());
    a_deviations[0] = a_f(a_a0.upper(), a_a1.nominal(), a_a2.nominal()) - nominal;
    a_deviations[1] = a_f(a_a0.nominal(), a_a1.upper(), a_a2.nominal()) - nominal;
    a_deviations[2] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.upper()) - nominal;
    return nominal;
  }

  template<typename F, typename T,
           typename N0, typename U0,
           typename N1, typename U1,
           typename N2, typename U2,
           typename N3, typename U3>
  static auto _propagate_error(F                        a_f,
                               static_vector<T, 4>&     a_deviations,
                               const uncertain<N0, U0>& a_a0,
                               const uncertain<N1, U1>& a_a1,
                               const uncertain<N2, U2>& a_a2,
                               const uncertain<N3, U3>& a_a3)
  {
    auto nominal    = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal());
    a_deviations[0] = a_f(a_a0.upper(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal()) - nominal;
    a_deviations[1] = a_f(a_a0.nominal(), a_a1.upper(), a_a2.nominal(), a_a3.nominal()) - nominal;
    a_deviations[2] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.upper(), a_a3.nominal()) - nominal;
    a_deviations[3] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.upper()) - nominal;
    return nominal;
  }

  template<typename F, typename T,
           typename N0, typename U0,
           typename N1, typename U1,
           typename N2, typename U2,
           typename N3, typename U3,
           typename N4, typename U4>
  static auto _propagate_error(F                        a_f,
                               static_vector<T, 5>&     a_deviations,
                               const uncertain<N0, U0>& a_a0,
                               const uncertain<N1, U1>& a_a1,
                               const uncertain<N2, U2>& a_a2,
                               const uncertain<N3, U3>& a_a3,
                               const uncertain<N4, U4>& a_a4)
  {
    auto nominal    = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal());
    a_deviations[0] = a_f(a_a0.upper(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal()) - nominal;
    a_deviations[1] = a_f(a_a0.nominal(), a_a1.upper(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal()) - nominal;
    a_deviations[2] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.upper(), a_a3.nominal(), a_a4.nominal()) - nominal;
    a_deviations[3] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.upper(), a_a4.nominal()) - nominal;
    a_deviations[4] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.upper()) - nominal;
    return nominal;
  }

  template<typename F, typename T,
           typename N0, typename U0,
           typename N1, typename U1,
           typename N2, typename U2,
           typename N3, typename U3,
           typename N4, typename U4,
           typename N5, typename U5>
  static auto _propagate_error(F                        a_f,
                               static_vector<T, 6>&     a_deviations,
                               const uncertain<N0, U0>& a_a0,
                               const uncertain<N1, U1>& a_a1,
                               const uncertain<N2, U2>& a_a2,
                               const uncertain<N3, U3>& a_a3,
                               const uncertain<N4, U4>& a_a4,
                               const uncertain<N5, U5>& a_a5)
  {
    auto nominal    = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal());
    a_deviations[0] = a_f(a_a0.upper(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal()) - nominal;
    a_deviations[1] = a_f(a_a0.nominal(), a_a1.upper(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal()) - nominal;
    a_deviations[2] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.upper(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal()) - nominal;
    a_deviations[3] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.upper(), a_a4.nominal(), a_a5.nominal()) - nominal;
    a_deviations[4] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.upper(), a_a5.nominal()) - nominal;
    a_deviations[5] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.upper()) - nominal;
    return nominal;
  }

  template<typename F, typename T,
           typename N0, typename U0,
           typename N1, typename U1,
           typename N2, typename U2,
           typename N3, typename U3,
           typename N4, typename U4,
           typename N5, typename U5,
           typename N6, typename U6>
  static auto _propagate_error(F                        a_f,
                               static_vector<T, 7>&     a_deviations,
                               const uncertain<N0, U0>& a_a0,
                               const uncertain<N1, U1>& a_a1,
                               const uncertain<N2, U2>& a_a2,
                               const uncertain<N3, U3>& a_a3,
                               const uncertain<N4, U4>& a_a4,
                               const uncertain<N5, U5>& a_a5,
                               const uncertain<N6, U6>& a_a6)
  {
    auto nominal    = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal());
    a_deviations[0] = a_f(a_a0.upper(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal()) - nominal;
    a_deviations[1] = a_f(a_a0.nominal(), a_a1.upper(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal()) - nominal;
    a_deviations[2] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.upper(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal()) - nominal;
    a_deviations[3] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.upper(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal()) - nominal;
    a_deviations[4] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.upper(), a_a5.nominal(), a_a6.nominal()) - nominal;
    a_deviations[5] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.upper(), a_a6.nominal()) - nominal;
    a_deviations[6] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.upper()) - nominal;
    return nominal;
  }

  template<typename F, typename T,
           typename N0, typename U0,
           typename N1, typename U1,
           typename N2, typename U2,
           typename N3, typename U3,
           typename N4, typename U4,
           typename N5, typename U5,
           typename N6, typename U6,
           typename N7, typename U7>
  static auto _propagate_error(F                        a_f,
                               static_vector<T, 8>&     a_deviations,
                               const uncertain<N0, U0>& a_a0,
                               const uncertain<N1, U1>& a_a1,
                               const uncertain<N2, U2>& a_a2,
                               const uncertain<N3, U3>& a_a3,
                               const uncertain<N4, U4>& a_a4,
                               const uncertain<N5, U5>& a_a5,
                               const uncertain<N6, U6>& a_a6,
                               const uncertain<N7, U7>& a_a7)
  {
    auto nominal    = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal());
    a_deviations[0] = a_f(a_a0.upper(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal()) - nominal;
    a_deviations[1] = a_f(a_a0.nominal(), a_a1.upper(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal()) - nominal;
    a_deviations[2] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.upper(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal()) - nominal;
    a_deviations[3] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.upper(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal()) - nominal;
    a_deviations[4] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.upper(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal()) - nominal;
    a_deviations[5] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.upper(), a_a6.nominal(), a_a7.nominal()) - nominal;
    a_deviations[6] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.upper(), a_a7.nominal()) - nominal;
    a_deviations[7] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.upper()) - nominal;
    return nominal;
  }

  template<typename F, typename T,
           typename N0, typename U0,
           typename N1, typename U1,
           typename N2, typename U2,
           typename N3, typename U3,
           typename N4, typename U4,
           typename N5, typename U5,
           typename N6, typename U6,
           typename N7, typename U7,
           typename N8, typename U8>
  static auto _propagate_error(F                        a_f,
                               static_vector<T, 9>&     a_deviations,
                               const uncertain<N0, U0>& a_a0,
                               const uncertain<N1, U1>& a_a1,
                               const uncertain<N2, U2>& a_a2,
                               const uncertain<N3, U3>& a_a3,
                               const uncertain<N4, U4>& a_a4,
                               const uncertain<N5, U5>& a_a5,
                               const uncertain<N6, U6>& a_a6,
                               const uncertain<N7, U7>& a_a7,
                               const uncertain<N8, U8>& a_a8)
  {
    auto nominal    = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal());
    a_deviations[0] = a_f(a_a0.upper(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal()) - nominal;
    a_deviations[1] = a_f(a_a0.nominal(), a_a1.upper(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal()) - nominal;
    a_deviations[2] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.upper(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal()) - nominal;
    a_deviations[3] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.upper(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal()) - nominal;
    a_deviations[4] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.upper(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal()) - nominal;
    a_deviations[5] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.upper(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal()) - nominal;
    a_deviations[6] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.upper(), a_a7.nominal(), a_a8.nominal()) - nominal;
    a_deviations[7] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.upper(), a_a8.nominal()) - nominal;
    a_deviations[8] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.upper()) - nominal;
    return nominal;
  }

  template<typename F, typename T,
           typename N0, typename U0,
           typename N1, typename U1,
           typename N2, typename U2,
           typename N3, typename U3,
           typename N4, typename U4,
           typename N5, typename U5,
           typename N6, typename U6,
           typename N7, typename U7,
           typename N8, typename U8,
           typename N9, typename U9>
  static auto _propagate_error(F                        a_f,
                               static_vector<T, 10>&    a_deviations,
                               const uncertain<N0, U0>& a_a0,
                               const uncertain<N1, U1>& a_a1,
                               const uncertain<N2, U2>& a_a2,
                               const uncertain<N3, U3>& a_a3,
                               const uncertain<N4, U4>& a_a4,
                               const uncertain<N5, U5>& a_a5,
                               const uncertain<N6, U6>& a_a6,
                               const uncertain<N7, U7>& a_a7,
                               const uncertain<N8, U8>& a_a8,
                               const uncertain<N9, U9>& a_a9)
  {
    auto nominal    = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal());
    a_deviations[0] = a_f(a_a0.upper(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal()) - nominal;
    a_deviations[1] = a_f(a_a0.nominal(), a_a1.upper(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal()) - nominal;
    a_deviations[2] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.upper(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal()) - nominal;
    a_deviations[3] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.upper(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal()) - nominal;
    a_deviations[4] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.upper(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal()) - nominal;
    a_deviations[5] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.upper(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal()) - nominal;
    a_deviations[6] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.upper(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal()) - nominal;
    a_deviations[7] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.upper(), a_a8.nominal(), a_a9.nominal()) - nominal;
    a_deviations[8] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.upper(), a_a9.nominal()) - nominal;
    a_deviations[9] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.upper()) - nominal;
    return nominal;
  }

  template<typename F, typename T,
           typename N0, typename U0,
           typename N1, typename U1,
           typename N2, typename U2,
           typename N3, typename U3,
           typename N4, typename U4,
           typename N5, typename U5,
           typename N6, typename U6,
           typename N7, typename U7,
           typename N8, typename U8,
           typename N9, typename U9,
           typename N10, typename U10>
  static auto _propagate_error(F                          a_f,
                               static_vector<T, 11>&      a_deviations,
                               const uncertain<N0, U0>&   a_a0,
                               const uncertain<N1, U1>&   a_a1,
                               const uncertain<N2, U2>&   a_a2,
                               const uncertain<N3, U3>&   a_a3,
                               const uncertain<N4, U4>&   a_a4,
                               const uncertain<N5, U5>&   a_a5,
                               const uncertain<N6, U6>&   a_a6,
                               const uncertain<N7, U7>&   a_a7,
                               const uncertain<N8, U8>&   a_a8,
                               const uncertain<N9, U9>&   a_a9,
                               const uncertain<N10, U10>& a_a10)
  {
    auto nominal     = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal());
    a_deviations[0]  = a_f(a_a0.upper(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal()) - nominal;
    a_deviations[1]  = a_f(a_a0.nominal(), a_a1.upper(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal()) - nominal;
    a_deviations[2]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.upper(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal()) - nominal;
    a_deviations[3]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.upper(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal()) - nominal;
    a_deviations[4]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.upper(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal()) - nominal;
    a_deviations[5]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.upper(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal()) - nominal;
    a_deviations[6]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.upper(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal()) - nominal;
    a_deviations[7]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.upper(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal()) - nominal;
    a_deviations[8]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.upper(), a_a9.nominal(), a_a10.nominal()) - nominal;
    a_deviations[9]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.upper(), a_a10.nominal()) - nominal;
    a_deviations[10] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.upper()) - nominal;
    return nominal;
  }

  template<typename F, typename T,
           typename N0, typename U0,
           typename N1, typename U1,
           typename N2, typename U2,
           typename N3, typename U3,
           typename N4, typename U4,
           typename N5, typename U5,
           typename N6, typename U6,
           typename N7, typename U7,
           typename N8, typename U8,
           typename N9, typename U9,
           typename N10, typename U10,
           typename N11, typename U11>
  static auto _propagate_error(F                          a_f,
                               static_vector<T, 12>&      a_deviations,
                               const uncertain<N0, U0>&   a_a0,
                               const uncertain<N1, U1>&   a_a1,
                               const uncertain<N2, U2>&   a_a2,
                               const uncertain<N3, U3>&   a_a3,
                               const uncertain<N4, U4>&   a_a4,
                               const uncertain<N5, U5>&   a_a5,
                               const uncertain<N6, U6>&   a_a6,
                               const uncertain<N7, U7>&   a_a7,
                               const uncertain<N8, U8>&   a_a8,
                               const uncertain<N9, U9>&   a_a9,
                               const uncertain<N10, U10>& a_a10,
                               const uncertain<N11, U11>& a_a11)
  {
    auto nominal     = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal());
    a_deviations[0]  = a_f(a_a0.upper(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal()) - nominal;
    a_deviations[1]  = a_f(a_a0.nominal(), a_a1.upper(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal()) - nominal;
    a_deviations[2]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.upper(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal()) - nominal;
    a_deviations[3]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.upper(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal()) - nominal;
    a_deviations[4]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.upper(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal()) - nominal;
    a_deviations[5]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.upper(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal()) - nominal;
    a_deviations[6]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.upper(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal()) - nominal;
    a_deviations[7]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.upper(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal()) - nominal;
    a_deviations[8]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.upper(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal()) - nominal;
    a_deviations[9]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.upper(), a_a10.nominal(), a_a11.nominal()) - nominal;
    a_deviations[10] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.upper(), a_a11.nominal()) - nominal;
    a_deviations[11] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.upper()) - nominal;
    return nominal;
  }

  template<typename F, typename T,
           typename N0, typename U0,
           typename N1, typename U1,
           typename N2, typename U2,
           typename N3, typename U3,
           typename N4, typename U4,
           typename N5, typename U5,
           typename N6, typename U6,
           typename N7, typename U7,
           typename N8, typename U8,
           typename N9, typename U9,
           typename N10, typename U10,
           typename N11, typename U11,
           typename N12, typename U12>
  static auto _propagate_error(F                          a_f,
                               static_vector<T, 13>&      a_deviations,
                               const uncertain<N0, U0>&   a_a0,
                               const uncertain<N1, U1>&   a_a1,
                               const uncertain<N2, U2>&   a_a2,
                               const uncertain<N3, U3>&   a_a3,
                               const uncertain<N4, U4>&   a_a4,
                               const uncertain<N5, U5>&   a_a5,
                               const uncertain<N6, U6>&   a_a6,
                               const uncertain<N7, U7>&   a_a7,
                               const uncertain<N8, U8>&   a_a8,
                               const uncertain<N9, U9>&   a_a9,
                               const uncertain<N10, U10>& a_a10,
                               const uncertain<N11, U11>& a_a11,
                               const uncertain<N12, U12>& a_a12)
  {
    auto nominal     = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal());
    a_deviations[0]  = a_f(a_a0.upper(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal()) - nominal;
    a_deviations[1]  = a_f(a_a0.nominal(), a_a1.upper(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal()) - nominal;
    a_deviations[2]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.upper(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal()) - nominal;
    a_deviations[3]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.upper(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal()) - nominal;
    a_deviations[4]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.upper(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal()) - nominal;
    a_deviations[5]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.upper(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal()) - nominal;
    a_deviations[6]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.upper(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal()) - nominal;
    a_deviations[7]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.upper(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal()) - nominal;
    a_deviations[8]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.upper(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal()) - nominal;
    a_deviations[9]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.upper(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal()) - nominal;
    a_deviations[10] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.upper(), a_a11.nominal(), a_a12.nominal()) - nominal;
    a_deviations[11] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.upper(), a_a12.nominal()) - nominal;
    a_deviations[12] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.upper()) - nominal;
    return nominal;
  }

  template<typename F, typename T,
           typename N0, typename U0,
           typename N1, typename U1,
           typename N2, typename U2,
           typename N3, typename U3,
           typename N4, typename U4,
           typename N5, typename U5,
           typename N6, typename U6,
           typename N7, typename U7,
           typename N8, typename U8,
           typename N9, typename U9,
           typename N10, typename U10,
           typename N11, typename U11,
           typename N12, typename U12,
           typename N13, typename U13>
  static auto _propagate_error(F                          a_f,
                               static_vector<T, 14>&      a_deviations,
                               const uncertain<N0, U0>&   a_a0,
                               const uncertain<N1, U1>&   a_a1,
                               const uncertain<N2, U2>&   a_a2,
                               const uncertain<N3, U3>&   a_a3,
                               const uncertain<N4, U4>&   a_a4,
                               const uncertain<N5, U5>&   a_a5,
                               const uncertain<N6, U6>&   a_a6,
                               const uncertain<N7, U7>&   a_a7,
                               const uncertain<N8, U8>&   a_a8,
                               const uncertain<N9, U9>&   a_a9,
                               const uncertain<N10, U10>& a_a10,
                               const uncertain<N11, U11>& a_a11,
                               const uncertain<N12, U12>& a_a12,
                               const uncertain<N13, U13>& a_a13)
  {
    auto nominal     = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal());
    a_deviations[0]  = a_f(a_a0.upper(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal()) - nominal;
    a_deviations[1]  = a_f(a_a0.nominal(), a_a1.upper(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal()) - nominal;
    a_deviations[2]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.upper(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal()) - nominal;
    a_deviations[3]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.upper(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal()) - nominal;
    a_deviations[4]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.upper(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal()) - nominal;
    a_deviations[5]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.upper(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal()) - nominal;
    a_deviations[6]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.upper(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal()) - nominal;
    a_deviations[7]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.upper(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal()) - nominal;
    a_deviations[8]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.upper(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal()) - nominal;
    a_deviations[9]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.upper(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal()) - nominal;
    a_deviations[10] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.upper(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal()) - nominal;
    a_deviations[11] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.upper(), a_a12.nominal(), a_a13.nominal()) - nominal;
    a_deviations[12] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.upper(), a_a13.nominal()) - nominal;
    a_deviations[13] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.upper()) - nominal;
    return nominal;
  }

  template<typename F, typename T,
           typename N0, typename U0,
           typename N1, typename U1,
           typename N2, typename U2,
           typename N3, typename U3,
           typename N4, typename U4,
           typename N5, typename U5,
           typename N6, typename U6,
           typename N7, typename U7,
           typename N8, typename U8,
           typename N9, typename U9,
           typename N10, typename U10,
           typename N11, typename U11,
           typename N12, typename U12,
           typename N13, typename U13,
           typename N14, typename U14>
  static auto _propagate_error(F                          a_f,
                               static_vector<T, 15>&      a_deviations,
                               const uncertain<N0, U0>&   a_a0,
                               const uncertain<N1, U1>&   a_a1,
                               const uncertain<N2, U2>&   a_a2,
                               const uncertain<N3, U3>&   a_a3,
                               const uncertain<N4, U4>&   a_a4,
                               const uncertain<N5, U5>&   a_a5,
                               const uncertain<N6, U6>&   a_a6,
                               const uncertain<N7, U7>&   a_a7,
                               const uncertain<N8, U8>&   a_a8,
                               const uncertain<N9, U9>&   a_a9,
                               const uncertain<N10, U10>& a_a10,
                               const uncertain<N11, U11>& a_a11,
                               const uncertain<N12, U12>& a_a12,
                               const uncertain<N13, U13>& a_a13,
                               const uncertain<N14, U14>& a_a14)
  {
    auto nominal     = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal());
    a_deviations[0]  = a_f(a_a0.upper(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal()) - nominal;
    a_deviations[1]  = a_f(a_a0.nominal(), a_a1.upper(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal()) - nominal;
    a_deviations[2]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.upper(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal()) - nominal;
    a_deviations[3]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.upper(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal()) - nominal;
    a_deviations[4]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.upper(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal()) - nominal;
    a_deviations[5]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.upper(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal()) - nominal;
    a_deviations[6]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.upper(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal()) - nominal;
    a_deviations[7]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.upper(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal()) - nominal;
    a_deviations[8]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.upper(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal()) - nominal;
    a_deviations[9]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.upper(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal()) - nominal;
    a_deviations[10] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.upper(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal()) - nominal;
    a_deviations[11] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.upper(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal()) - nominal;
    a_deviations[12] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.upper(), a_a13.nominal(), a_a14.nominal()) - nominal;
    a_deviations[13] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.upper(), a_a14.nominal()) - nominal;
    a_deviations[14] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.upper()) - nominal;
    return nominal;
  }

  template<typename F, typename T,
           typename N0, typename U0,
           typename N1, typename U1,
           typename N2, typename U2,
           typename N3, typename U3,
           typename N4, typename U4,
           typename N5, typename U5,
           typename N6, typename U6,
           typename N7, typename U7,
           typename N8, typename U8,
           typename N9, typename U9,
           typename N10, typename U10,
           typename N11, typename U11,
           typename N12, typename U12,
           typename N13, typename U13,
           typename N14, typename U14,
           typename N15, typename U15>
  static auto _propagate_error(F                          a_f,
                               static_vector<T, 16>&      a_deviations,
                               const uncertain<N0, U0>&   a_a0,
                               const uncertain<N1, U1>&   a_a1,
                               const uncertain<N2, U2>&   a_a2,
                               const uncertain<N3, U3>&   a_a3,
                               const uncertain<N4, U4>&   a_a4,
                               const uncertain<N5, U5>&   a_a5,
                               const uncertain<N6, U6>&   a_a6,
                               const uncertain<N7, U7>&   a_a7,
                               const uncertain<N8, U8>&   a_a8,
                               const uncertain<N9, U9>&   a_a9,
                               const uncertain<N10, U10>& a_a10,
                               const uncertain<N11, U11>& a_a11,
                               const uncertain<N12, U12>& a_a12,
                               const uncertain<N13, U13>& a_a13,
                               const uncertain<N14, U14>& a_a14,
                               const uncertain<N15, U15>& a_a15)
  {
    auto nominal     = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal());
    a_deviations[0]  = a_f(a_a0.upper(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal()) - nominal;
    a_deviations[1]  = a_f(a_a0.nominal(), a_a1.upper(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal()) - nominal;
    a_deviations[2]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.upper(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal()) - nominal;
    a_deviations[3]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.upper(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal()) - nominal;
    a_deviations[4]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.upper(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal()) - nominal;
    a_deviations[5]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.upper(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal()) - nominal;
    a_deviations[6]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.upper(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal()) - nominal;
    a_deviations[7]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.upper(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal()) - nominal;
    a_deviations[8]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.upper(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal()) - nominal;
    a_deviations[9]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.upper(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal()) - nominal;
    a_deviations[10] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.upper(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal()) - nominal;
    a_deviations[11] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.upper(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal()) - nominal;
    a_deviations[12] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.upper(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal()) - nominal;
    a_deviations[13] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.upper(), a_a14.nominal(), a_a15.nominal()) - nominal;
    a_deviations[14] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.upper(), a_a15.nominal()) - nominal;
    a_deviations[15] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.upper()) - nominal;
    return nominal;
  }

  template<typename F, typename T,
           typename N0, typename U0,
           typename N1, typename U1,
           typename N2, typename U2,
           typename N3, typename U3,
           typename N4, typename U4,
           typename N5, typename U5,
           typename N6, typename U6,
           typename N7, typename U7,
           typename N8, typename U8,
           typename N9, typename U9,
           typename N10, typename U10,
           typename N11, typename U11,
           typename N12, typename U12,
           typename N13, typename U13,
           typename N14, typename U14,
           typename N15, typename U15,
           typename N16, typename U16>
  static auto _propagate_error(F                          a_f,
                               static_vector<T, 17>&      a_deviations,
                               const uncertain<N0, U0>&   a_a0,
                               const uncertain<N1, U1>&   a_a1,
                               const uncertain<N2, U2>&   a_a2,
                               const uncertain<N3, U3>&   a_a3,
                               const uncertain<N4, U4>&   a_a4,
                               const uncertain<N5, U5>&   a_a5,
                               const uncertain<N6, U6>&   a_a6,
                               const uncertain<N7, U7>&   a_a7,
                               const uncertain<N8, U8>&   a_a8,
                               const uncertain<N9, U9>&   a_a9,
                               const uncertain<N10, U10>& a_a10,
                               const uncertain<N11, U11>& a_a11,
                               const uncertain<N12, U12>& a_a12,
                               const uncertain<N13, U13>& a_a13,
                               const uncertain<N14, U14>& a_a14,
                               const uncertain<N15, U15>& a_a15,
                               const uncertain<N16, U16>& a_a16)
  {
    auto nominal     = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal(), a_a16.nominal());
    a_deviations[0]  = a_f(a_a0.upper(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal(), a_a16.nominal()) - nominal;
    a_deviations[1]  = a_f(a_a0.nominal(), a_a1.upper(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal(), a_a16.nominal()) - nominal;
    a_deviations[2]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.upper(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal(), a_a16.nominal()) - nominal;
    a_deviations[3]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.upper(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal(), a_a16.nominal()) - nominal;
    a_deviations[4]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.upper(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal(), a_a16.nominal()) - nominal;
    a_deviations[5]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.upper(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal(), a_a16.nominal()) - nominal;
    a_deviations[6]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.upper(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal(), a_a16.nominal()) - nominal;
    a_deviations[7]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.upper(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal(), a_a16.nominal()) - nominal;
    a_deviations[8]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.upper(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal(), a_a16.nominal()) - nominal;
    a_deviations[9]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.upper(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal(), a_a16.nominal()) - nominal;
    a_deviations[10] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.upper(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal(), a_a16.nominal()) - nominal;
    a_deviations[11] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.upper(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal(), a_a16.nominal()) - nominal;
    a_deviations[12] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.upper(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal(), a_a16.nominal()) - nominal;
    a_deviations[13] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.upper(), a_a14.nominal(), a_a15.nominal(), a_a16.nominal()) - nominal;
    a_deviations[14] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.upper(), a_a15.nominal(), a_a16.nominal()) - nominal;
    a_deviations[15] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.upper(), a_a16.nominal()) - nominal;
    a_deviations[16] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal(), a_a16.upper()) - nominal;
    return nominal;
  }

  template<typename F, typename T,
           typename N0, typename U0,
           typename N1, typename U1,
           typename N2, typename U2,
           typename N3, typename U3,
           typename N4, typename U4,
           typename N5, typename U5,
           typename N6, typename U6,
           typename N7, typename U7,
           typename N8, typename U8,
           typename N9, typename U9,
           typename N10, typename U10,
           typename N11, typename U11,
           typename N12, typename U12,
           typename N13, typename U13,
           typename N14, typename U14,
           typename N15, typename U15,
           typename N16, typename U16,
           typename N17, typename U17>
  static auto _propagate_error(F                          a_f,
                               static_vector<T, 18>&      a_deviations,
                               const uncertain<N0, U0>&   a_a0,
                               const uncertain<N1, U1>&   a_a1,
                               const uncertain<N2, U2>&   a_a2,
                               const uncertain<N3, U3>&   a_a3,
                               const uncertain<N4, U4>&   a_a4,
                               const uncertain<N5, U5>&   a_a5,
                               const uncertain<N6, U6>&   a_a6,
                               const uncertain<N7, U7>&   a_a7,
                               const uncertain<N8, U8>&   a_a8,
                               const uncertain<N9, U9>&   a_a9,
                               const uncertain<N10, U10>& a_a10,
                               const uncertain<N11, U11>& a_a11,
                               const uncertain<N12, U12>& a_a12,
                               const uncertain<N13, U13>& a_a13,
                               const uncertain<N14, U14>& a_a14,
                               const uncertain<N15, U15>& a_a15,
                               const uncertain<N16, U16>& a_a16,
                               const uncertain<N17, U17>& a_a17)
  {
    auto nominal     = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal(), a_a16.nominal(), a_a17.nominal());
    a_deviations[0]  = a_f(a_a0.upper(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal(), a_a16.nominal(), a_a17.nominal()) - nominal;
    a_deviations[1]  = a_f(a_a0.nominal(), a_a1.upper(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal(), a_a16.nominal(), a_a17.nominal()) - nominal;
    a_deviations[2]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.upper(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal(), a_a16.nominal(), a_a17.nominal()) - nominal;
    a_deviations[3]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.upper(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal(), a_a16.nominal(), a_a17.nominal()) - nominal;
    a_deviations[4]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.upper(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal(), a_a16.nominal(), a_a17.nominal()) - nominal;
    a_deviations[5]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.upper(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal(), a_a16.nominal(), a_a17.nominal()) - nominal;
    a_deviations[6]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.upper(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal(), a_a16.nominal(), a_a17.nominal()) - nominal;
    a_deviations[7]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.upper(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal(), a_a16.nominal(), a_a17.nominal()) - nominal;
    a_deviations[8]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.upper(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal(), a_a16.nominal(), a_a17.nominal()) - nominal;
    a_deviations[9]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.upper(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal(), a_a16.nominal(), a_a17.nominal()) - nominal;
    a_deviations[10] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.upper(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal(), a_a16.nominal(), a_a17.nominal()) - nominal;
    a_deviations[11] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.upper(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal(), a_a16.nominal(), a_a17.nominal()) - nominal;
    a_deviations[12] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.upper(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal(), a_a16.nominal(), a_a17.nominal()) - nominal;
    a_deviations[13] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.upper(), a_a14.nominal(), a_a15.nominal(), a_a16.nominal(), a_a17.nominal()) - nominal;
    a_deviations[14] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.upper(), a_a15.nominal(), a_a16.nominal(), a_a17.nominal()) - nominal;
    a_deviations[15] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.upper(), a_a16.nominal(), a_a17.nominal()) - nominal;
    a_deviations[16] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal(), a_a16.upper(), a_a17.nominal()) - nominal;
    a_deviations[17] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal(), a_a16.nominal(), a_a17.upper()) - nominal;
    return nominal;
  }

  template<typename F, typename T,
           typename N0, typename U0,
           typename N1, typename U1,
           typename N2, typename U2,
           typename N3, typename U3,
           typename N4, typename U4,
           typename N5, typename U5,
           typename N6, typename U6,
           typename N7, typename U7,
           typename N8, typename U8,
           typename N9, typename U9,
           typename N10, typename U10,
           typename N11, typename U11,
           typename N12, typename U12,
           typename N13, typename U13,
           typename N14, typename U14,
           typename N15, typename U15,
           typename N16, typename U16,
           typename N17, typename U17,
           typename N18, typename U18>
  static auto _propagate_error(F                          a_f,
                               static_vector<T, 19>&      a_deviations,
                               const uncertain<N0, U0>&   a_a0,
                               const uncertain<N1, U1>&   a_a1,
                               const uncertain<N2, U2>&   a_a2,
                               const uncertain<N3, U3>&   a_a3,
                               const uncertain<N4, U4>&   a_a4,
                               const uncertain<N5, U5>&   a_a5,
                               const uncertain<N6, U6>&   a_a6,
                               const uncertain<N7, U7>&   a_a7,
                               const uncertain<N8, U8>&   a_a8,
                               const uncertain<N9, U9>&   a_a9,
                               const uncertain<N10, U10>& a_a10,
                               const uncertain<N11, U11>& a_a11,
                               const uncertain<N12, U12>& a_a12,
                               const uncertain<N13, U13>& a_a13,
                               const uncertain<N14, U14>& a_a14,
                               const uncertain<N15, U15>& a_a15,
                               const uncertain<N16, U16>& a_a16,
                               const uncertain<N17, U17>& a_a17,
                               const uncertain<N18, U18>& a_a18)
  {
    auto nominal     = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal(), a_a16.nominal(), a_a17.nominal(), a_a18.nominal());
    a_deviations[0]  = a_f(a_a0.upper(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal(), a_a16.nominal(), a_a17.nominal(), a_a18.nominal()) - nominal;
    a_deviations[1]  = a_f(a_a0.nominal(), a_a1.upper(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal(), a_a16.nominal(), a_a17.nominal(), a_a18.nominal()) - nominal;
    a_deviations[2]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.upper(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal(), a_a16.nominal(), a_a17.nominal(), a_a18.nominal()) - nominal;
    a_deviations[3]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.upper(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal(), a_a16.nominal(), a_a17.nominal(), a_a18.nominal()) - nominal;
    a_deviations[4]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.upper(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal(), a_a16.nominal(), a_a17.nominal(), a_a18.nominal()) - nominal;
    a_deviations[5]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.upper(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal(), a_a16.nominal(), a_a17.nominal(), a_a18.nominal()) - nominal;
    a_deviations[6]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.upper(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal(), a_a16.nominal(), a_a17.nominal(), a_a18.nominal()) - nominal;
    a_deviations[7]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.upper(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal(), a_a16.nominal(), a_a17.nominal(), a_a18.nominal()) - nominal;
    a_deviations[8]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.upper(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal(), a_a16.nominal(), a_a17.nominal(), a_a18.nominal()) - nominal;
    a_deviations[9]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.upper(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal(), a_a16.nominal(), a_a17.nominal(), a_a18.nominal()) - nominal;
    a_deviations[10] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.upper(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal(), a_a16.nominal(), a_a17.nominal(), a_a18.nominal()) - nominal;
    a_deviations[11] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.upper(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal(), a_a16.nominal(), a_a17.nominal(), a_a18.nominal()) - nominal;
    a_deviations[12] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.upper(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal(), a_a16.nominal(), a_a17.nominal(), a_a18.nominal()) - nominal;
    a_deviations[13] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.upper(), a_a14.nominal(), a_a15.nominal(), a_a16.nominal(), a_a17.nominal(), a_a18.nominal()) - nominal;
    a_deviations[14] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.upper(), a_a15.nominal(), a_a16.nominal(), a_a17.nominal(), a_a18.nominal()) - nominal;
    a_deviations[15] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.upper(), a_a16.nominal(), a_a17.nominal(), a_a18.nominal()) - nominal;
    a_deviations[16] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal(), a_a16.upper(), a_a17.nominal(), a_a18.nominal()) - nominal;
    a_deviations[17] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal(), a_a16.nominal(), a_a17.upper(), a_a18.nominal()) - nominal;
    a_deviations[18] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal(), a_a16.nominal(), a_a17.nominal(), a_a18.upper()) - nominal;
    return nominal;
  }

  template<typename F, typename T,
           typename N0, typename U0,
           typename N1, typename U1,
           typename N2, typename U2,
           typename N3, typename U3,
           typename N4, typename U4,
           typename N5, typename U5,
           typename N6, typename U6,
           typename N7, typename U7,
           typename N8, typename U8,
           typename N9, typename U9,
           typename N10, typename U10,
           typename N11, typename U11,
           typename N12, typename U12,
           typename N13, typename U13,
           typename N14, typename U14,
           typename N15, typename U15,
           typename N16, typename U16,
           typename N17, typename U17,
           typename N18, typename U18,
           typename N19, typename U19>
  static auto _propagate_error(F                          a_f,
                               static_vector<T, 20>&      a_deviations,
                               const uncertain<N0, U0>&   a_a0,
                               const uncertain<N1, U1>&   a_a1,
                               const uncertain<N2, U2>&   a_a2,
                               const uncertain<N3, U3>&   a_a3,
                               const uncertain<N4, U4>&   a_a4,
                               const uncertain<N5, U5>&   a_a5,
                               const uncertain<N6, U6>&   a_a6,
                               const uncertain<N7, U7>&   a_a7,
                               const uncertain<N8, U8>&   a_a8,
                               const uncertain<N9, U9>&   a_a9,
                               const uncertain<N10, U10>& a_a10,
                               const uncertain<N11, U11>& a_a11,
                               const uncertain<N12, U12>& a_a12,
                               const uncertain<N13, U13>& a_a13,
                               const uncertain<N14, U14>& a_a14,
                               const uncertain<N15, U15>& a_a15,
                               const uncertain<N16, U16>& a_a16,
                               const uncertain<N17, U17>& a_a17,
                               const uncertain<N18, U18>& a_a18,
                               const uncertain<N19, U19>& a_a19)
  {
    auto nominal     = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal(), a_a16.nominal(), a_a17.nominal(), a_a18.nominal(), a_a19.nominal());
    a_deviations[0]  = a_f(a_a0.upper(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal(), a_a16.nominal(), a_a17.nominal(), a_a18.nominal(), a_a19.nominal()) - nominal;
    a_deviations[1]  = a_f(a_a0.nominal(), a_a1.upper(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal(), a_a16.nominal(), a_a17.nominal(), a_a18.nominal(), a_a19.nominal()) - nominal;
    a_deviations[2]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.upper(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal(), a_a16.nominal(), a_a17.nominal(), a_a18.nominal(), a_a19.nominal()) - nominal;
    a_deviations[3]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.upper(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal(), a_a16.nominal(), a_a17.nominal(), a_a18.nominal(), a_a19.nominal()) - nominal;
    a_deviations[4]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.upper(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal(), a_a16.nominal(), a_a17.nominal(), a_a18.nominal(), a_a19.nominal()) - nominal;
    a_deviations[5]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.upper(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal(), a_a16.nominal(), a_a17.nominal(), a_a18.nominal(), a_a19.nominal()) - nominal;
    a_deviations[6]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.upper(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal(), a_a16.nominal(), a_a17.nominal(), a_a18.nominal(), a_a19.nominal()) - nominal;
    a_deviations[7]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.upper(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal(), a_a16.nominal(), a_a17.nominal(), a_a18.nominal(), a_a19.nominal()) - nominal;
    a_deviations[8]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.upper(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal(), a_a16.nominal(), a_a17.nominal(), a_a18.nominal(), a_a19.nominal()) - nominal;
    a_deviations[9]  = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.upper(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal(), a_a16.nominal(), a_a17.nominal(), a_a18.nominal(), a_a19.nominal()) - nominal;
    a_deviations[10] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.upper(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal(), a_a16.nominal(), a_a17.nominal(), a_a18.nominal(), a_a19.nominal()) - nominal;
    a_deviations[11] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.upper(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal(), a_a16.nominal(), a_a17.nominal(), a_a18.nominal(), a_a19.nominal()) - nominal;
    a_deviations[12] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.upper(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal(), a_a16.nominal(), a_a17.nominal(), a_a18.nominal(), a_a19.nominal()) - nominal;
    a_deviations[13] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.upper(), a_a14.nominal(), a_a15.nominal(), a_a16.nominal(), a_a17.nominal(), a_a18.nominal(), a_a19.nominal()) - nominal;
    a_deviations[14] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.upper(), a_a15.nominal(), a_a16.nominal(), a_a17.nominal(), a_a18.nominal(), a_a19.nominal()) - nominal;
    a_deviations[15] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.upper(), a_a16.nominal(), a_a17.nominal(), a_a18.nominal(), a_a19.nominal()) - nominal;
    a_deviations[16] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal(), a_a16.upper(), a_a17.nominal(), a_a18.nominal(), a_a19.nominal()) - nominal;
    a_deviations[17] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal(), a_a16.nominal(), a_a17.upper(), a_a18.nominal(), a_a19.nominal()) - nominal;
    a_deviations[18] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal(), a_a16.nominal(), a_a17.nominal(), a_a18.upper(), a_a19.nominal()) - nominal;
    a_deviations[19] = a_f(a_a0.nominal(), a_a1.nominal(), a_a2.nominal(), a_a3.nominal(), a_a4.nominal(), a_a5.nominal(), a_a6.nominal(), a_a7.nominal(), a_a8.nominal(), a_a9.nominal(), a_a10.nominal(), a_a11.nominal(), a_a12.nominal(), a_a13.nominal(), a_a14.nominal(), a_a15.nominal(), a_a16.nominal(), a_a17.nominal(), a_a18.nominal(), a_a19.upper()) - nominal;
    return nominal;
  }

  // END GENERATED CODE
};

}  // namespace libUncertainty
