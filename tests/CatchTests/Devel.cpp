#define CATCH_CONFIG_ENABLE_BENCHMARKING
#include "catch.hpp"

#include <iostream>
#include <sstream>

#include <BoostUnitDefinitions/Units.hpp>
#include <boost/any.hpp>
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

#include <map>
#include <variant>

using namespace boost::units;
using namespace libUncertainty;

TEST_CASE("Usage")
{
  SECTION("Sigfig rounding")
  {
    SECTION("doubles")
    {
      CHECK(sigfig_round(1.23456, 1) == Approx(1));
      CHECK(sigfig_round(1.23456, 2) == Approx(1.2));
      CHECK(sigfig_round(1.23456, 3) == Approx(1.23));
      CHECK(sigfig_round(1.23456, 4) == Approx(1.235));
    }
    SECTION("boost quantities")
    {
      quantity<t::m> x = 1.23456 * i::m;

      CHECK(sigfig_round(x, 1).value() == Approx(1));
      CHECK(sigfig_round(x, 2).value() == Approx(1.2));
      CHECK(sigfig_round(x, 3).value() == Approx(1.23));
      CHECK(sigfig_round(x, 4).value() == Approx(1.235));
    }
  }

  SECTION("Constructing uncertain<...> variables")
  {
    SECTION("Default")
    {
      uncertain<double> x;

      x.nominal(2.1);
      x.uncertainty(0.1);

      CHECK(x.nominal() == Approx(2.1));
      CHECK(x.uncertainty() == Approx(0.1));
      CHECK(x.upper() == Approx(2.2));
      CHECK(x.lower() == Approx(2.0));

      CHECK(x.relative_uncertainty() == Approx(0.1 / 2.1));
    }
    SECTION("With nominal value and uncertainty")
    {
      uncertain<int> x(2, 1);
      CHECK(x.nominal() == 2);
      CHECK(x.uncertainty() == 1);
    }
    SECTION("With nominal value only")
    {
      uncertain<int> x(2);
      CHECK(x.nominal() == 2);
      CHECK(x.uncertainty() == 0);
    }
    SECTION("With nominal initializer list")
    {
      SECTION("to ctor")
      {
        uncertain<int> x{2, 1};
        CHECK(x.nominal() == 2);
        CHECK(x.uncertainty() == 1);
      }
      SECTION("assignment")
      {
        uncertain<int> x = {2, 1};
        CHECK(x.nominal() == 2);
        CHECK(x.uncertainty() == 1);
      }
    }
    SECTION("With array")
    {
      std::array<int, 2> data = {2, 1};
      uncertain<int>     x(data);
      CHECK(x.nominal() == 2);
      CHECK(x.uncertainty() == 1);
    }
    SECTION("With vector")
    {
      std::vector<int> data = {2, 1};
      uncertain<int>   x(data);
      CHECK(x.nominal() == 2);
      CHECK(x.uncertainty() == 1);
    }
  }
  SECTION("init constructed")
  {
    uncertain<double> x(2.1, 0.1);

    CHECK(x.nominal() == Approx(2.1));
    CHECK(x.uncertainty() == Approx(0.1));
    CHECK(x.upper() == Approx(2.2));
    CHECK(x.lower() == Approx(2.0));
  }

  SECTION("iostream output")
  {
    uncertain<int>    x(10, 2);
    std::stringstream out;

    out << x;

    CHECK(out.str() == "10 +/- 2");
  }

  SECTION("Error Propagation")
  {
    SECTION("passing lvalue")
    {
      uncertain<double> x(2, 0.2);
      auto              y = basic_error_propagator::propagate_error([](double x) { return x * x; }, x);

      CHECK(y.nominal() == Approx(4));
      CHECK(y.uncertainty() == Approx(2.2 * 2.2 - 4));

      y = basic_error_propagator::propagate_error([](double x) { return -x * x; }, x);

      CHECK(y.nominal() == Approx(-4));
      CHECK(y.uncertainty() == Approx(2.2 * 2.2 - 4));
    }
    SECTION("passing rvalue")
    {
      auto y = basic_error_propagator::propagate_error([](double x) { return x * x; }, uncertain<double>(2, 0.2));

      CHECK(y.nominal() == Approx(4));
      CHECK(y.uncertainty() == Approx(2.2 * 2.2 - 4));
    }
  }

  SECTION("Boost.Units Integration")
  {
    SECTION("Length quantity")
    {
      uncertain<quantity<t::cm>> L(2 * i::cm, 0.1 * i::cm);
      uncertain<quantity<t::cm>> W(4 * i::cm, 0.2 * i::cm);

      auto A = basic_error_propagator::propagate_error([](quantity<t::cm> L, quantity<t::cm> W) { return L * W; }, L, W);

      CHECK(A.nominal().value() == Approx(8));
      CHECK(A.uncertainty().value() == Approx(0.5656854));
      CHECK(quantity<t::m_p2>(A.nominal()).value() == Approx(0.0008));
      CHECK(quantity<t::mm_p2>(A.uncertainty()).value() == Approx(56.56854));
    }
    SECTION("Different units")
    {
      SECTION("with full type")
      {
        uncertain<quantity<t::cm>, quantity<t::mm>> L(2 * i::cm, 1 * i::mm);
        uncertain<quantity<t::cm>, quantity<t::mm>> W(4 * i::cm, 2 * i::mm);

        auto A = basic_error_propagator::propagate_error([](quantity<t::cm> L, quantity<t::cm> W) { return L * W; }, L, W);

        CHECK(A.nominal().value() == Approx(8));
        CHECK(A.uncertainty().value() == Approx(0.5656854));
        CHECK(quantity<t::m_p2>(A.nominal()).value() == Approx(0.0008));
        CHECK(quantity<t::mm_p2>(A.uncertainty()).value() == Approx(56.56854));
      }

      SECTION("with helper function")
      {
        auto L = make_uncertain(2 * i::cm, 1 * i::mm);
        auto W = make_uncertain(4 * i::cm, 2 * i::mm);

        auto A = basic_error_propagator::propagate_error([](quantity<t::cm> L, quantity<t::cm> W) { return L * W; }, L, W);

        CHECK(A.nominal().value() == Approx(8));
        CHECK(A.uncertainty().value() == Approx(0.5656854));
        CHECK(quantity<t::m_p2>(A.nominal()).value() == Approx(0.0008));
        CHECK(quantity<t::mm_p2>(A.uncertainty()).value() == Approx(56.56854));
      }
    }
  }

  SECTION("Normalization")
  {
    SECTION("doubles")
    {
      auto              x = make_uncertain(1.23456, 0.98765);
      std::stringstream out;

      out << x.normalize();
      CHECK(out.str() == "1 +/- 1");
      out.str("");
      out << x.normalize(1);
      CHECK(out.str() == "1 +/- 1");
      out.str("");
      out << x.normalize(2);
      CHECK(out.str() == "1.23 +/- 0.99");
      out.str("");
      out << x.normalize(3);
      CHECK(out.str() == "1.235 +/- 0.988");
      out.str("");
      out << x.normalize(4);
      CHECK(out.str() == "1.2346 +/- 0.9877");
    }
    SECTION("Boost quantities")
    {
      SECTION("Same unit for nominal and uncertainty")
      {
        auto              x = make_uncertain(1.23456 * i::m, 0.98765 * i::m);
        std::stringstream out;

        out << x.normalize();
        CHECK(out.str() == "1 m +/- 1 m");
        out.str("");
        out << x.normalize(1);
        CHECK(out.str() == "1 m +/- 1 m");
        out.str("");
        out << x.normalize(2);
        CHECK(out.str() == "1.23 m +/- 0.99 m");
        out.str("");
        out << x.normalize(3);
        CHECK(out.str() == "1.235 m +/- 0.988 m");
        out.str("");
        out << x.normalize(4);
        CHECK(out.str() == "1.2346 m +/- 0.9877 m");
      }

      SECTION("Different unit for nominal and uncertainty")
      {
        auto              x = make_uncertain(1.23456 * i::m, 98.765 * i::cm);
        std::stringstream out;

        out << x.normalize();
        CHECK(out.str() == "1 m +/- 100 cm");
        out.str("");
        out << x.normalize(1);
        CHECK(out.str() == "1 m +/- 100 cm");
        out.str("");
        out << x.normalize(2);
        CHECK(out.str() == "1.23 m +/- 99 cm");
        out.str("");
        out << x.normalize(3);
        CHECK(out.str() == "1.235 m +/- 98.8 cm");
        out.str("");
        out << x.normalize(4);
        CHECK(out.str() == "1.2346 m +/- 98.77 cm");
      }
    }
  }
}

TEST_CASE("Correlations")
{
  SECTION("Correlation matrix")
  {
    correlation_matrix<double> mat(3);
    mat(0, 1) = 0.1;
    mat(0, 2) = 0.2;
    mat(1, 2) = 0.3;

    CHECK(mat(1, 0) == Approx(0.1));
  }

  SECTION("Error propagation w/ correlation matrix.")
  {
    boost::numeric::ublas::matrix<double> corr(2, 2);
    corr(0, 0) = 1;
    corr(1, 0) = 1;
    corr(0, 1) = 1;
    corr(1, 1) = 1;

    uncertain<double> x(2, 0.1), y(3, 0.1);

    auto z = basic_error_propagator::propagate_error([](double a, double b) { return b - a; }, corr, x, y);

    // a and b are directly correlated, so there is no uncertainy in their difference.
    // the nominal and upper values should be the same.
    CHECK(z.nominal() == Approx(1));
    CHECK(z.upper() == Approx(1));
  }

  SECTION("add_correlation_coefficients Mixin")
  {
    add_correlation_coefficient_array<uncertain<double>> x({2.2, 0.1});
    add_correlation_coefficient_array<uncertain<double>> y = {3.3, 0.2};

    CHECK(x.get_correlation_coefficients().size() == 0);
    x.set_correlation_coefficient_array_size(3);
    x.get_correlation_coefficients()[0] = 1;
    x.get_correlation_coefficients()[1] = 2;
    x.get_correlation_coefficients()[2] = 3;

    auto z = x;

    CHECK(x.get_correlation_coefficients().size() == 3);
    CHECK(x.get_correlation_coefficients()[0] == Approx(1));
    CHECK(x.get_correlation_coefficients()[1] == Approx(2));
    CHECK(x.get_correlation_coefficients()[2] == Approx(3));

    CHECK(z.get_correlation_coefficients().size() == 3);
    CHECK(z.get_correlation_coefficients()[0] == Approx(1));
    CHECK(z.get_correlation_coefficients()[1] == Approx(2));
    CHECK(z.get_correlation_coefficients()[2] == Approx(3));
  }

  SECTION("Correlation store")
  {
    add_id<uncertain<double>> x,y,z;
    correlation_store<double> store;

    store.add(x, y, 0.1);
    store.add(y, z, 0.2);

    CHECK_THROWS(store.add(x, y, 0.1));

    CHECK(store.get(x, y) == Approx(0.1));
    CHECK(store.get(y, z) == Approx(0.2));

    auto& global_store = get_global_correlation_store();

    global_store.add(x, y, 0.1);
    global_store.add(x, z, 0.2);

    CHECK_THROWS(global_store.add(y, x, 0.1));

    CHECK(global_store.get(x, y) == Approx(0.1));
    CHECK(global_store.get(x, z) == Approx(0.2));

    global_store.set(x,y,0.5);
    global_store.set(x,z,-0.5);
    global_store.set(z,y,1);

    CHECK(global_store.get(y,x) == Approx(0.5) );
    CHECK(global_store.get(z,x) == Approx(-0.5) );
    CHECK(global_store.get(z,y) == Approx(1) );

  }

  SECTION("Error propagation w/ correlation")
  {
    SECTION("Doubles")
    {
      uncertain<double> x(1, 0.1), y(2, 0.2);

      SECTION("Independent")
      {
        auto z_and_correlation = basic_error_propagator::propagate_error_and_correlation([](double a, double b) { return a + b; }, x, y);

        double unc = sqrt(0.1 * 0.1 + 0.2 * 0.2);
        CHECK(z_and_correlation.nominal() == Approx(3));
        CHECK(z_and_correlation.upper() == Approx(3 + unc));
        CHECK(z_and_correlation.get_correlation_coefficients().size() == 2);
        CHECK(z_and_correlation.get_correlation_coefficient(0) == Approx(0.1 / unc));
        CHECK(z_and_correlation.get_correlation_coefficient(1) == Approx(0.2 / unc));
      }

      SECTION("with Correlation Matrix")
      {
        correlation_matrix<double> corr(2);
        corr(0, 1) = -1;

        auto z_and_correlation = basic_error_propagator::propagate_error_and_correlation([](double a, double b) { return a + b; }, corr, x, y);

        double unc = 0.1;
        CHECK(z_and_correlation.nominal() == Approx(3));
        CHECK(z_and_correlation.upper() == Approx(3 + unc));
        CHECK(z_and_correlation.get_correlation_coefficients().size() == 2);
        CHECK(z_and_correlation.get_correlation_coefficients().size() == 2);
        CHECK(z_and_correlation.get_correlation_coefficient(0) == Approx((0.1 + 0.2 * -1) / unc));
        CHECK(z_and_correlation.get_correlation_coefficient(1) == Approx((0.2 + 0.1 * -1) / unc));
      }
      SECTION("with Correlation Store")
      {
        correlation_store<double> cstore;
        add_id<uncertain<double>> xx;
        add_id<uncertain<double>> yy;
        xx = x;
        yy = y;

        cstore.add(xx, yy, -1);

        auto z = basic_error_propagator::propagate_error([](double a, double b) { return a + b; }, cstore, xx, yy);

        double unc = 0.1;
        CHECK(z.nominal() == Approx(3));
        CHECK(z.upper() == Approx(3 + unc));
        CHECK(cstore.get(z, xx) == Approx((0.1 + 0.2 * -1) / unc));
        CHECK(cstore.get(z, yy) == Approx((0.2 + 0.1 * -1) / unc));

      }

    }

    SECTION("boost quantities")
    {
      uncertain<quantity<t::cm>> x(1 * i::cm, 0.1 * i::cm), y(2 * i::cm, 0.2 * i::cm);
      SECTION("Independent")
      {
        SECTION("Addition")
        {
          auto z_and_correlation = basic_error_propagator::propagate_error_and_correlation([](auto a, auto b) { return a + b; }, x, y);

          double unc = sqrt(0.1 * 0.1 + 0.2 * 0.2);
          CHECK(z_and_correlation.nominal().value() == Approx(3));
          CHECK(z_and_correlation.upper().value() == Approx(3 + unc));
          CHECK(z_and_correlation.get_correlation_coefficients().size() == 2);
          CHECK(z_and_correlation.get_correlation_coefficients().size() == 2);
          CHECK(z_and_correlation.get_correlation_coefficient(0) == Approx(0.1 / unc));
          CHECK(z_and_correlation.get_correlation_coefficient(1) == Approx(0.2 / unc));
        }

        SECTION("Multiplyication")
        {
          auto z_and_correlation = basic_error_propagator::propagate_error_and_correlation([](auto a, auto b) { return a * b; }, x, y);

          double unc = sqrt(0.2 * 0.2 + 0.2 * 0.2);
          CHECK(z_and_correlation.nominal().value() == Approx(2));
          CHECK(z_and_correlation.upper().value() == Approx(2 + unc));
          CHECK(z_and_correlation.get_correlation_coefficients().size() == 2);
          CHECK(z_and_correlation.get_correlation_coefficient(0) == Approx(0.2 / unc));
          CHECK(z_and_correlation.get_correlation_coefficient(1) == Approx(0.2 / unc));
        }
      }
      SECTION("with Correlation Matrix")
      {
        correlation_matrix<double> corr(2);
        corr(0, 1) = -1;

        auto z_and_correlation = basic_error_propagator::propagate_error_and_correlation([](auto a, auto b) { return a + b; }, corr, x, y);

        double unc = 0.1;
        CHECK(z_and_correlation.nominal().value() == Approx(3));
        CHECK(z_and_correlation.upper().value() == Approx(3 + unc));
        CHECK(z_and_correlation.get_correlation_coefficients().size() == 2);
        CHECK(z_and_correlation.get_correlation_coefficients().size() == 2);
        CHECK(z_and_correlation.get_correlation_coefficient(0) == Approx((0.1 + 0.2 * -1) / unc));
        CHECK(z_and_correlation.get_correlation_coefficient(1) == Approx((0.2 + 0.1 * -1) / unc));
      }
    }
  }
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

TEST_CASE("id'ed variables")
{
  SECTION("getting unique ids")
  {
    CHECK(get_uniq_id() > 0);
    CHECK(get_uniq_id() > 1);
  }

  SECTION("id'ed uncertain double")
  {
    add_id<uncertain<double>> x;
    x = uncertain<double>(2, 0.2);
    add_id<uncertain<float>>  y(1, 0.1);
    add_id<uncertain<double>> z;
    z = x;

    CHECK(y.get_id() > x.get_id());
    CHECK(z.get_id() == x.get_id());

    CHECK(z.nominal() == Approx(2));
    CHECK(z.upper() == Approx(2.2));

    z.new_id();
    CHECK(z.get_id() > y.get_id());

    CHECK( get_id(x) == x.get_id() );
    CHECK( get_id(1) == 0 );
    double u;
    uncertain<double> v;
    CHECK( get_id(u) == 0 );
    CHECK( get_id(v) == 0 );

  }
}
