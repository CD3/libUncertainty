#define CATCH_CONFIG_ENABLE_BENCHMARKING
#include "catch.hpp"

#include <iostream>
#include <sstream>
#include <libUncertainty/uncertain.hpp>
#include <BoostUnitDefinitions/Units.hpp>

using namespace boost::units;
using namespace libUncertainty;

TEST_CASE("uncertain<...> tests")
{
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
  SECTION("make_uncertain with data")
  {
    SECTION("doubles")
    {
      std::vector<double> vals{0.431, 0.603, 0.504, 0.581, 0.588, 0.644, 0.595, 0.534, 0.563, 0.578};

      auto t = make_uncertain(vals.begin(), vals.end());
      CHECK( t.nominal() == Approx(0.5620999999999999) );
      CHECK( t.uncertainty() == Approx(0.059718692402146764/sqrt(10)) );

      t = make_uncertain(vals.begin(), vals.end(), tags::use_stdev_for_error{});
      CHECK( t.nominal() == Approx(0.5620999999999999) );
      CHECK( t.uncertainty() == Approx(0.059718692402146764) );
    }

  }
}
