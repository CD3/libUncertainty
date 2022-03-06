#include "catch.hpp"

#include <iostream>
#include <sstream>

#include <libUncertainty/propagate.hpp>
#include <libUncertainty/uncertain.hpp>
#include <libUncertainty/utils.hpp>

#include <BoostUnitDefinitions/Units.hpp>
/**
 * This file is used for developement. As new classes are created, small tests
 * are written here so that we can try to compile and use them.
 */

using namespace boost::units;
using namespace libUncertainty;

TEST_CASE("General", "[devel]")
{
  SECTION("Sigfig rounding")
  {
    CHECK(sigfig_round(1.23456, 1) == Approx(1));
    CHECK(sigfig_round(1.23456, 2) == Approx(1.2));
    CHECK(sigfig_round(1.23456, 3) == Approx(1.23));
    CHECK(sigfig_round(1.23456, 4) == Approx(1.235));
  }

  SECTION("Constructing")
  {
    SECTION("Default")
    {
      uncertain<double> x;
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

  SECTION("Usage")
  {
    SECTION("default constructed")
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
    SECTION("init constructed")
    {
      uncertain<double> x(2.1, 0.1);

      CHECK(x.nominal() == Approx(2.1));
      CHECK(x.uncertainty() == Approx(0.1));
      CHECK(x.upper() == Approx(2.2));
      CHECK(x.lower() == Approx(2.0));
    }

    SECTION("With boost::quantity")
    {
      // uncertain<quantity<t::m>> x;
    }
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
      auto              y = error_propagators::basic::propagate_error([](double x) { return x * x; }, x);

      CHECK(y.nominal() == Approx(4));
      CHECK(y.uncertainty() == Approx(2.2 * 2.2 - 4));

      y = error_propagators::basic::propagate_error([](double x) { return -x * x; }, x);

      CHECK(y.nominal() == Approx(-4));
      CHECK(y.uncertainty() == Approx(2.2 * 2.2 - 4));
    }
    SECTION("passing rvalue")
    {
      auto y = error_propagators::basic::propagate_error([](double x) { return x * x; }, uncertain<double>(2, 0.2));

      CHECK(y.nominal() == Approx(4));
      CHECK(y.uncertainty() == Approx(2.2 * 2.2 - 4));
    }
  }

  SECTION("Boost.Units Integration")
  {
    SECTION("Length quantity")
    {
      uncertain<quantity<t::cm>> L(2*i::cm,0.1*i::cm);
      uncertain<quantity<t::cm>> W(4*i::cm,0.2*i::cm);

      auto A = error_propagators::basic::propagate_error([](quantity<t::cm> L, quantity<t::cm> W) { return L * W; }, L, W);

      CHECK( A.nominal().value() == Approx(8) );
      CHECK( A.uncertainty().value() == Approx(0.5656854) );
      CHECK( quantity<t::m_p2>(A.nominal()).value() == Approx(0.0008) );
      CHECK( quantity<t::mm_p2>(A.uncertainty()).value() == Approx(56.56854) );


    }
    SECTION("Different units")
    {
      SECTION("with full type")
      {
      uncertain<quantity<t::cm>,quantity<t::mm>> L(2*i::cm,1*i::mm);
      uncertain<quantity<t::cm>,quantity<t::mm>> W(4*i::cm,2*i::mm);

      auto A = error_propagators::basic::propagate_error([](quantity<t::cm> L, quantity<t::cm> W) { return L * W; }, L, W);

      CHECK( A.nominal().value() == Approx(8) );
      CHECK( A.uncertainty().value() == Approx(0.5656854) );
      CHECK( quantity<t::m_p2>(A.nominal()).value() == Approx(0.0008) );
      CHECK( quantity<t::mm_p2>(A.uncertainty()).value() == Approx(56.56854) );
      }

      SECTION("with helper function")
      {
      auto L = make_uncertain(2*i::cm,1*i::mm);
      auto W = make_uncertain(4*i::cm,2*i::mm);

      auto A = error_propagators::basic::propagate_error([](quantity<t::cm> L, quantity<t::cm> W) { return L * W; }, L, W);

      CHECK( A.nominal().value() == Approx(8) );
      CHECK( A.uncertainty().value() == Approx(0.5656854) );
      CHECK( quantity<t::m_p2>(A.nominal()).value() == Approx(0.0008) );
      CHECK( quantity<t::mm_p2>(A.uncertainty()).value() == Approx(56.56854) );
      }
    }
  }
}
