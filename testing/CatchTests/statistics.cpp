#define CATCH_CONFIG_ENABLE_BENCHMARKING
#include <iostream>

#include <BoostUnitDefinitions/Units.hpp>

#include <catch2/catch_all.hpp>
#include <libUncertainty/statistics.hpp>

using namespace libUncertainty;
using namespace boost::units;
using namespace Catch;

TEST_CASE("Builtin statistic calculation tests")
{
  SECTION("statistic calculations")
  {
    SECTION("doubles")
    {
      std::vector<double> vals{0.431, 0.603, 0.504, 0.581, 0.588, 0.644, 0.595, 0.534, 0.563, 0.578};

      double avg = average(vals.begin(), vals.end());
      CHECK(avg == Approx(0.5620999999999999));

      double var = variance(vals.begin(), vals.end());
      CHECK(var == Approx(0.0035663222222222218));

      var = variance(vals.begin(), vals.end(), 0);
      CHECK(var == Approx(0.003209689999999999));

      double std = standard_deviation(vals.begin(), vals.end());
      CHECK(std == Approx(0.059718692402146764));

      std = standard_deviation(vals.begin(), vals.end(), 0);
      CHECK(std == Approx(0.05665412606333275));

      double err = standard_error_of_the_mean(vals.begin(), vals.end());
      CHECK(err == Approx(0.059718692402146764 / sqrt(10)));
    }
    SECTION("quantities")
    {
      std::vector<quantity<t::s>> vals{0.431 * i::s, 0.603 * i::s, 0.504 * i::s, 0.581 * i::s, 0.588 * i::s, 0.644 * i::s, 0.595 * i::s, 0.534 * i::s, 0.563 * i::s, 0.578 * i::s};

      quantity<t::s> avg = average(vals.begin(), vals.end());
      CHECK(avg.value() == Approx(0.5620999999999999));
      using s_p2 = typename decltype(avg * avg)::unit_type;

      quantity<s_p2> var = variance(vals.begin(), vals.end());
      CHECK(var.value() == Approx(0.0035663222222222218));

      var = variance(vals.begin(), vals.end(), 0);
      CHECK(var.value() == Approx(0.003209689999999999));

      quantity<t::s> std = standard_deviation(vals.begin(), vals.end());
      CHECK(std.value() == Approx(0.059718692402146764));
      CHECK(quantity<t::ms>(std).value() == Approx(59.718692402146764));

      std = standard_deviation(vals.begin(), vals.end(), 0);
      CHECK(std.value() == Approx(0.05665412606333275));
      CHECK(quantity<t::ms>(std).value() == Approx(56.65412606333275));
    }
  }
}
