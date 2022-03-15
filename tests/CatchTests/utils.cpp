
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

  SECTION("Uncertainy free functions")
  {
    double x = 10;
    uncertain<double> y(20,2);

    CHECK( !is_uncertain(x) );
    CHECK(  is_uncertain(y) );
    CHECK( !is_uncertain(1) );
    CHECK(  is_uncertain(uncertain<double>(3,0.1)) );

    CHECK( get_nominal(x) == Approx(10));
    CHECK( get_uncertainty(x) == Approx(0).scale(1));
    CHECK( get_upper(x) == Approx(10));
    CHECK( get_lower(x) == Approx(10));

    CHECK( get_nominal(y) == Approx(20));
    CHECK( get_uncertainty(y) == Approx(2));
    CHECK( get_upper(y) == Approx(22));
    CHECK( get_lower(y) == Approx(18));

    CHECK( get_nominal(1) == Approx(1));
    CHECK( get_uncertainty(1) == Approx(0).scale(1));
    CHECK( get_upper(1) == Approx(1));
    CHECK( get_lower(1) == Approx(1));

    CHECK( get_nominal(uncertain<double>(3,0.1)) == Approx(3));
    CHECK( get_uncertainty(uncertain<double>(3,0.1)) == Approx(0.1));
    CHECK( get_upper(uncertain<double>(3,0.1)) == Approx(3.1));
    CHECK( get_lower(uncertain<double>(3,0.1)) == Approx(2.9));



  }
}
