
#define CATCH_CONFIG_ENABLE_BENCHMARKING
#include "catch.hpp"

#include <libUncertainty/utils.hpp>
#include <libUncertainty/uncertain.hpp>
#include <BoostUnitDefinitions/Units.hpp>


using namespace boost::units;
using namespace libUncertainty;

TEST_CASE("Utility functions/classes Tests")
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

  SECTION("getting unique ids")
  {
    auto id1 = get_uniq_id();
    auto id2 = get_uniq_id();
    auto id3 = get_uniq_id();

    CHECK(id1 > 0);
    CHECK(id2 != id1);
    CHECK(id3 != id1);
    CHECK(id3 != id2);
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

    CHECK(get_id(x) == x.get_id());
    CHECK(get_id(1) == 0);
    double            u;
    uncertain<double> v;
    CHECK(get_id(u) == 0);
    CHECK(get_id(v) == 0);
  }
}
