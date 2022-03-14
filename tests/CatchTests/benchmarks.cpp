#define CATCH_CONFIG_ENABLE_BENCHMARKING
#include "catch.hpp"

#include <BoostUnitDefinitions/Units.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/type_traits/function_traits.hpp>

#include <libUncertainty/correlation.hpp>
#include <libUncertainty/propagate.hpp>
#include <libUncertainty/uncertain.hpp>
#include <libUncertainty/utils.hpp>
// clang-format off
#include <uncertainties/ureal.hpp>
#include <uncertainties/impl.hpp>
#include <uncertainties/io.hpp>
#include <uncertainties/math.hpp>
// clang-format on
/**
 * This file is used for developement. As new classes are created, small tests
 * are written here so that we can try to compile and use them.
 */

using namespace boost::units;
using namespace libUncertainty;

TEST_CASE("Memory Usage")
{
  CHECK( sizeof( uncertain<double> ) == 2*sizeof(double) );
  CHECK( sizeof( add_id<uncertain<double>> ) == 2*sizeof(double)+sizeof(size_t) );
  CHECK( sizeof( uncertainties::udouble ) > sizeof(uncertain<double>));
  CHECK( sizeof( uncertainties::udouble ) > sizeof(add_id<uncertain<double>>));
}

TEST_CASE("Bencharmks", "[.][benchmarks]")
{
  SECTION("Error Propagation")
  {
    auto              f = [](double x, double y, double z) { return sin(x) * cos(y) * tan(z); };
    uncertain<double> x(M_PI / 2, 0.01), y(M_PI, 0.01), z(M_PI / 4, 0.01);
    BENCHMARK("Normal Call")
    {
      return f(x.nominal(), y.nominal(), z.nominal());
    };
    BENCHMARK("Error Propagation w/o Correlations")
    {
      return basic_error_propagator::propagate_error(f, x, y, z);
    };

    auto u = basic_error_propagator::propagate_error(f, x, y, z);
    CHECK(u.nominal() == Approx(-1));
    CHECK(u.uncertainty() == Approx(0.0202028));
  }

  SECTION("uncertainties-cpp comparison")
  {
    SECTION("uncertainties-cpp calculations")
    {
      uncertainties::udouble alpha(0.200, 0.003), gamma(0.100, 0.001);
      BENCHMARK("Solid angle")
      {
        return 2 * M_PI * (1 - uncertainties::cos(alpha / 2));
      };
      uncertainties::udouble Omega = 2 * M_PI * (1 - uncertainties::cos(alpha / 2));
      std::stringstream      out;
      out << Omega;

      CHECK(out.str() == "(3.14 ± 0.09)e-2");
    }
    SECTION("libUncertainty calculations")
    {
      uncertain<double> alpha(.200, 0.003), gamma(0.100, 0.001);
      BENCHMARK("Solid angle")
      {
        return basic_error_propagator::propagate_error([](double alpha) { return 2 * M_PI * (1 - cos(alpha / 2)); }, alpha);
      };
      auto              Omega = basic_error_propagator::propagate_error([](double alpha) { return 2 * M_PI * (1 - cos(alpha / 2)); }, alpha);
      std::stringstream out;
      out << Omega;
      CHECK(out.str() == "0.0313898 +/- 0.000947941");
      out.str("");
      out << Omega.normalize();
      CHECK(out.str() == "0.0314 +/- 0.0009");
    }
  }
}
