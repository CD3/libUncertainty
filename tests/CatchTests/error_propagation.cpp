#define CATCH_CONFIG_ENABLE_BENCHMARKING
#include "catch.hpp"

#include <BoostUnitDefinitions/Units.hpp>
#include <boost/any.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/type_traits/function_traits.hpp>

#include <libUncertainty/propagate.hpp>
#include <libUncertainty/uncertain.hpp>

using namespace boost::units;
using namespace libUncertainty;

TEST_CASE("Error Propagation Tests")
{
  SECTION("Basic Method")
  {
    SECTION("No correlation")
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
  }
  SECTION("Mixing uncertain with exact quantities")
  {
    SECTION("doubles")
    {
      double            x = 2;
      uncertain<double> y(2, 0.1);

      auto L = basic_error_propagator::propagate_error([](double a, double b) { return a + b; }, x, y);

      CHECK(L.nominal() == Approx(4));
      CHECK(L.uncertainty() == Approx(0.1));
    }
    SECTION("Lots of args")
    {
      auto L = basic_error_propagator::propagate_error([](double a, double b, double c, double d, double e, double f, double g) { return a + b + c + d + e + f + g; },
                                                       uncertain<double>{1, 0.1},
                                                       uncertain<double>{2, 0.2},
                                                       uncertain<double>{3, 0.3},
                                                       uncertain<double>{4, 0.4},
                                                       uncertain<double>{5, 0.5},
                                                       uncertain<double>{6, 0.6},
                                                       uncertain<double>{7, 0.7});

      CHECK(L.nominal() == Approx(28));
      CHECK(L.uncertainty() == Approx(1.1832159566));

      L = basic_error_propagator::propagate_error([](double a, double b, double c, double d, double e, double f, double g) { return a + b + c + d + e + f + g; },
                                                  1, uncertain<double>{2, 0.2}, 3, 4, 5, 6, 7);

      CHECK(L.nominal() == Approx(28));
      CHECK(L.uncertainty() == Approx(0.2));
    }

    SECTION("quantities")
    {
      auto h = make_uncertain(1.5*i::m, 1*i::cm);
      auto t = make_uncertain(0.562*i::s, 19*i::ms);

      auto f = [](quantity<t::m> h, quantity<t::s> t){ return 2*h/t/t; };

      auto g = basic_error_propagator::propagate_error(f, h, t);

      g = g.normalize();
      CHECK( g.nominal().value() == Approx(9.5) );
      CHECK( g.uncertainty().value() == Approx(0.6) );

    }
  }
}
