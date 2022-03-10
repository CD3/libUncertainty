#define CATCH_CONFIG_ENABLE_BENCHMARKING
#include "catch.hpp"

#include <iostream>
#include <sstream>

#include <BoostUnitDefinitions/Units.hpp>
#include <boost/any.hpp>
#include <boost/type_traits/function_traits.hpp>
#include <boost/numeric/ublas/matrix.hpp>

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

struct variable_store {
 private:
  std::map<const std::string, boost::any> m_storage;

 public:
  template<typename T>
  void add(const std::string& a_name)
  {
    m_storage[a_name] = T();
  }

  template<typename T>
  void add(const std::string& a_name, T a_t)
  {
    m_storage[a_name] = T(std::move(a_t));
  }

  template<typename T>
  T& get(const std::string& a_name)
  {
    return *boost::any_cast<T>(&m_storage.at(a_name));
  }
};

struct uncertainty_convention {
 private:
  variable_store m_store;

 public:
  variable_store&       get_variable_store() { return m_store; }
  const variable_store& get_variable_store() const { return m_store; }

  template<typename T>
  void add_variable(const std::string& a_name)
  {
    m_store.add<uncertain<T>>(a_name);
  }

  template<typename T>
  void add_variable(const std::string& a_name, T a_t)
  {
    m_store.add(a_name, uncertain<T>(std::move(a_t)));
  }

  template<typename T, typename U>
  void add_variable(const std::string& a_name, T a_t, U a_u)
  {
    m_store.add(a_name, uncertain<T, U>(std::move(a_t), std::move(a_u)));
  }

  template<typename T>
  uncertain<T>& get_variable(const std::string& a_name)
  {
    return m_store.get<uncertain<T>>(a_name);
  }

  template<typename T>
  const uncertain<T>& get_variable(const std::string& a_name) const
  {
    return m_store.get<uncertain<T>>(a_name);
  }

  template<typename FR, typename FA1, typename FA2>
  uncertain<FR> propagate_error(std::function<FR(FA1, FA2)> a_f, const std::string a_a1, const std::string a_a2)
  {
    auto nz  = a_f(this->get_variable<FA1>(a_a1).nominal(),
                  this->get_variable<FA2>(a_a2).nominal());
    auto dz1 = a_f(this->get_variable<FA1>(a_a1).upper(),
                   this->get_variable<FA2>(a_a2).nominal()) -
               nz;
    auto dz2 = a_f(this->get_variable<FA1>(a_a1).nominal(),
                   this->get_variable<FA2>(a_a2).upper()) -
               nz;
    auto dz = sqrt(dz1 * dz1 + dz2 * dz2);

    std::array<double, 2> correlation_coefficients;

    correlation_coefficients[0] = dz1 / dz;
    correlation_coefficients[1] = dz2 / dz;

    uncertain<FR> z(nz, dz);

    return {nz, dz};
  }
};

TEST_CASE("uncertainty_convention")
{
  SECTION("variable_store")
  {
    variable_store store;
    store.add<double>("x", 1.1);
    CHECK(store.get<double>("x") == Approx(1.1));
    store.get<double>("x") = 2.1;
    CHECK(store.get<double>("x") == Approx(2.1));
  }

  uncertainty_convention conv;
  conv.add_variable<double>("x", 2.1, 0.1);
  conv.add_variable<double>("y", 3.2, 0.2);

  {
    auto z = conv.propagate_error(std::function<double(double, double)>(std::plus<double>()), "x", "y");
    CHECK(z.nominal() == Approx(5.3));
    CHECK(z.uncertainty() == Approx(sqrt(0.05)));
  }

  {
    // auto z = conv.propagate_error([](double x, double y)->double{return x+y;}, "x", "y");
    // CHECK(z.nominal() == Approx(5.3));
  }
}

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
  boost::numeric::ublas::matrix<double> corr(2,2);
  corr(0,0) = 1;
  corr(1,0) = 1;
  corr(0,1) = 1;
  corr(1,1) = 1;

  uncertain<double> x(2,0.1),y(3,0.1);

  auto z = basic_error_propagator::propagate_error( [](double a, double b){ return b-a; }, corr, x, y);

  // a and b are directly correlated, so there is no uncertainy in their difference.
  // the nominal and upper values should be the same.
  CHECK( z.nominal() == Approx(1) );
  CHECK( z.upper() == Approx(1) );
}

TEST_CASE("Bencharmks")
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
      std::stringstream out;
      out << Omega;

      CHECK( out.str() == "(3.14 ± 0.09)e-2" );

    }
    SECTION("libUncertainty calculations")
    {
      uncertain<double> alpha(.200, 0.003), gamma(0.100, 0.001);
      BENCHMARK("Solid angle")
      {
        return basic_error_propagator::propagate_error([](double alpha) { return 2 * M_PI * (1 - cos(alpha / 2)); }, alpha);
      };
      auto Omega = basic_error_propagator::propagate_error([](double alpha) { return 2 * M_PI * (1 - cos(alpha / 2)); }, alpha);
      std::stringstream out;
      out << Omega;
      CHECK( out.str() == "0.0313898 +/- 0.000947941" );
      out.str("");
      out << Omega.normalize();
      CHECK( out.str() == "0.0314 +/- 0.0009" );

    }
  }
}
