#define CATCH_CONFIG_ENABLE_BENCHMARKING
#include "catch.hpp"

#include <BoostUnitDefinitions/Units.hpp>
#include <boost/any.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/type_traits/function_traits.hpp>

#include <libUncertainty/uncertain.hpp>
#include <libUncertainty/propagate.hpp>


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
}
