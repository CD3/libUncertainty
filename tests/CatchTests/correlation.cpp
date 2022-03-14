#define CATCH_CONFIG_ENABLE_BENCHMARKING
#include "catch.hpp"
#include <boost/numeric/ublas/matrix.hpp>

#include <libUncertainty/correlation.hpp>
#include <libUncertainty/uncertain.hpp>
#include <libUncertainty/propagate.hpp>
#include <BoostUnitDefinitions/Units.hpp>

using namespace boost::units;
using namespace libUncertainty;

TEST_CASE("Correlations Utilities")
{
  SECTION("Correlation matrix")
  {
    correlation_matrix<double> mat(3);
    CHECK(mat(0, 1) == Approx(0).scale(1));
    CHECK(mat(1, 0) == Approx(0).scale(1));
    CHECK(mat(0, 2) == Approx(0).scale(1));
    CHECK(mat(2, 0) == Approx(0).scale(1));
    CHECK(mat(0, 3) == Approx(0).scale(1));
    CHECK(mat(3, 0) == Approx(0).scale(1));
    CHECK(mat(0, 0) == Approx(1));
    CHECK(mat(1, 1) == Approx(1));
    CHECK(mat(2, 2) == Approx(1));

    mat(0, 1) = 0.1;
    mat(0, 2) = 0.2;
    mat(1, 2) = 0.3;

    CHECK(mat(0, 1) == Approx(0.1));
    CHECK(mat(1, 0) == Approx(0.1));
    CHECK(mat(0, 2) == Approx(0.2));
    CHECK(mat(2, 0) == Approx(0.2));
    CHECK(mat(0, 3) == Approx(0.3));
    CHECK(mat(3, 0) == Approx(0.3));
    CHECK(mat(0, 0) == Approx(1));
    CHECK(mat(1, 1) == Approx(1));
    CHECK(mat(2, 2) == Approx(1));
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
