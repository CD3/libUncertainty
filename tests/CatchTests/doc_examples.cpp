#define CATCH_CONFIG_ENABLE_BENCHMARKING
#include "catch.hpp"

#include <BoostUnitDefinitions/Units.hpp>

#include <libUncertainty/correlation.hpp>
#include <libUncertainty/propagate.hpp>
#include <libUncertainty/uncertain.hpp>
#include <libUncertainty/utils.hpp>
#include <libUncertainty/statistics.hpp>

using namespace libUncertainty;

double my_calculation(double x, double y)
{
  return 2 * x * y + 3 * x + 2 * y + 10;
}

TEST_CASE("README Examples")
{
  SECTION("The uncertain<...> class template")
  {
  }

  SECTION("Error Propgation without Correlation")
  {
    uncertain<double> x(1, 0.1);
    uncertain<double> y(3, 0.2);

    SECTION("With free function")
    {
      auto z = basic_error_propagator::propagate_error(my_calculation, x, y);
      z.nominal();      // 4
      z.uncertainty();  // 1.204159479

      CHECK(z.nominal() == Approx(25));
      CHECK(z.uncertainty() == Approx(1.204159479));
    }
    SECTION("With a lambda")
    {
      auto z = basic_error_propagator::propagate_error([](double x, double y) { return x + y; }, x, y);
      z.nominal();      // 4
      z.uncertainty();  // 0.2236067

      CHECK(z.nominal() == Approx(4));
      CHECK(z.uncertainty() == Approx(0.2236067));
    }
  }
  SECTION("Error Propgation with Correlation")
  {
    uncertain<double> x(4, 0.1);
    uncertain<double> y(3, 0.2);

    auto z = basic_error_propagator::propagate_error([](double a, double b) { return a + b; }, x, y);

    CHECK(z.nominal() == Approx(7));
    CHECK(z.uncertainty() == Approx(0.2236067977));

    correlation_matrix<double> corr(2);
    corr(0, 1) = -1;
    z          = basic_error_propagator::propagate_error([](double a, double b) { return a + b; }, corr, x, y);

    CHECK(z.nominal() == Approx(7));
    CHECK(z.uncertainty() == Approx(0.1));
  }

  SECTION("Error Propgation with Correlation using Store")
  {
    /*
    template<typename T>
    using tracked_uncertain<T> = add_id<uncertain<T>>;
    */

    correlation_store<double> store;

    add_id<uncertain<double>> x(4, 0.1);
    add_id<uncertain<double>> y(3, 0.2);

    store.set(x, y, -1);

    auto z = basic_error_propagator::propagate_error([](double a, double b) { return a + b; }, store, x, y);

    CHECK(z.nominal() == Approx(7));
    CHECK(z.uncertainty() == Approx(0.1));
    CHECK(store.get(z, x) == Approx(-1));
    CHECK(store.get(z, y) == Approx(1));

    auto& gstore = get_global_correlation_store();
  }

  SECTION("Gravity Example")
  {
    using namespace boost::units;
    std::vector<quantity<t::s>> time_data {
        0.431*i::s
      , 0.603*i::s
      , 0.504*i::s
      , 0.581*i::s
      , 0.588*i::s
      , 0.644*i::s
      , 0.595*i::s
      , 0.534*i::s
      , 0.563*i::s
      , 0.578*i::s
      };

    auto time_measurement = make_uncertain(time_data.begin(), time_data.end());
    auto height_measurement = make_uncertain( 1.5*i::m, 1*i::cm );

    auto calc_gravity = [](quantity<t::m> h, quantity<t::s> t) { return 2*h/t/t; };

    auto g = basic_error_propagator::propagate_error( calc_gravity, height_measurement, time_measurement );
    g = g.normalize();

    std::cout << "Measured g: " << g << std::endl;
    std::cout << "z-score: " << z_score( g, 9.81*i::m/i::s/i::s ) << std::endl;

    CHECK( g.nominal().value() == Approx(9.5) );
    CHECK( g.uncertainty().value() == Approx(0.6) );
  }
}
