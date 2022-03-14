# libUncertainty

A C++ library for working with uncertain quantities.

This is a _simple_ library, largely based on my Python library [pyErrorProp](https://github.com/CD3/pyErrorProp)
that provides some _basic_ utilities for working with uncertain quantities, including _basic_ error propagation.
It was written with [Boost.Units](https://www.boost.org/doc/libs/1_78_0/doc/html/boost_units.html) in mind, so
it supports uncertain quantities with units.

## Features

- Simple: the `uncertain<...>` class template is a thin wrapper around `std::pair<...>`.
- Piecemeal adoptions: you can just use the `uncertain<...>` class template to track uncertainty using your own methods, use the error propagation functions
  to propagate error with or with or without correlations, or use the
  correlation store to automatically track correlations between uncertain
  variables. Getting started is easy, and you don't have to commit to a new
  framework. Use it where it makes sense.
- Built with units support: the library is written and tested with support for the boost units library in mind. The library can propagate error through functions that take
  quantities (variables with units) as arguments.
- Non-intrusive: the library will propagate error through *arbitrary* functions, no special implementations or modification are necessary. Just create `uncertain<...>` variables and propagate
  error through your code.



## Installation

Build and install with CMake:
```
$ mkdir build
$ cd build
$ cmake ..
$ cmake --build .
$ cmake --build . --target install
```
Or, if you use conan
```
$ conan remote add cd3 https://cdc3.jfrog.io/artifactory/api/conan/default-conan
$ conan install libuncertainty/0.1@cd3/devel --build missing -g cmake_find_package
```


## Usage

### The `uncertain<...>` class template

An uncertain variable is created with the `uncertain<...>` class template. An uncertain variable
has a nominal value and an uncertainty.
```
#include <libUncertainty/uncertain.hpp>
...
using namespace libUncertainty;
...

uncertain<double> x(5.4,0.2);

x.nominal();     // 5.4
x.uncertainty(); // 0.2

std::cout << x << "\n"; // prints 5.4 +/- 0.2
```

Uncertain variables can use different types for the nominal value and uncertainty
```
uncertain<int,double> x(5,0.2);

x.nominal();     // 5
x.uncertainty(); // 0.2

std::cout << x << "\n"; // prints 5 +/- 0.2
```
This is useful when working with quantities with units. The nominal value and uncertainty can be specified with different units
```
#include <boost/units/quantity.hpp>
#include <BoostUnitDefinitions/Units.hpp> // include some convienient Boost.Unit definitions
using namespace boost::units;
...

uncertain<quantity<t::m>,quantity<t::cm>> L(12*i::m, 5*i::cm);

L.nominal();     // 12 m: a boost.units quantity with t::m as the unit
L.uncertainty(); // 5 cm: a boost.units quantity with t::cm as the unit
L.nominal().value();     // 12
L.uncertainty().value(); // 5

std::cout << L << "\n"; // prints 12 m +/- 5 cm
```

By itself, `uncertain<...>` does not provide any error propagation. It is just a simple container that provides some convenient methods
for working with uncertainty

```
uncertain<quantity<t::m>,quantity<t::cm>> L(1.2345*i::m, 2.34*i::cm);

L.upper();     // 1.2579
L.lower();     // 1.2111
L.normalize(); // returns a copy of the uncertain varaible with the
               // uncertainty rounded to one significant figure
               // and the nominal value rounded to the same decimal position.
               // 1.23 m +/- 2 cm

L.normalize(2); // same as above, but the uncertainy is rounded to two sigfigs.
               // 1.235 m +/- 2.3 cm

```
### Error Propagation (without Correlation)

Calculations with uncertain variables is not directly supported. Instead, error propagation (through arbitrary functions)
is done with a wrapper function.
```
#include <libUncertainty/propagate.hpp>
#include <libUncertainty/uncertain.hpp>
...
uncertain<double> x(1, 0.1);
uncertain<double> y(3, 0.2);

auto z = basic_error_propagator::propagate_error(my_calculation, x, y);
z.nominal();      // 4
z.uncertainty();  // 1.204159479
```
If you want to propagate error through a one-off calculation, you can use a lambda instead.
```
auto              z = basic_error_propagator::propagate_error([](double x, double y) { return x + y; }, x, y);
z.nominal();      // 4
z.uncertainty();  // 0.2236067
```
This makes the library *much* simpler than autograd-based libraries because we do not have to provide any special uncertain-aware functions.
It works with boost quantities too.
```
uncertain<quantity<t::cm>> L(2 * i::cm, 0.1 * i::cm);
uncertain<quantity<t::cm>> W(4 * i::cm, 0.2 * i::cm);

auto A = basic_error_propagator::propagate_error([](quantity<t::cm> L, quantity<t::cm> W) { return L * W; }, L, W);

A.nominal();      // 8 cm
A.uncertainty();  // 0.5656754 cm
```

### Error Propagation with Correlation

If the inputs to a function are correlated, the uncertainty in the result will depend on the correlation. You can tell the error propagator to
account for this by passing in a correlation matrix.

```
uncertain<double> x(4, 0.1);
uncertain<double> y(3, 0.2);

// uncorrelated inputs
auto z = basic_error_propagator::propagate_error([](double a, double b) { return a + b; }, x, y);

z.nominal();     //  7
z.uncertainty(); //  0.2236067977

correlation_matrix<double> corr(2);
corr(0, 1) = -1;

// the first and second inputs are corrleated
z = basic_error_propagator::propagate_error([](double a, double b) { return a + b; }, corr, x, y);

z.nominal();     // 7
z.uncertainty(); // 0.1
```

The correlation matrix is a square, symmetric matrix with whose elements give the correlation between arguments (the element (0,1) gives the correlation coefficient between the first and second arguments).
Since any variable is exactly correlated with itself, the diagonal elements of the correlation matrix will be one. The `correlation_matrix<...>` class template is a simple container for creating correlation matrices.
You only need to tell it the size, and then set the non-diagonal elements in the upper (or lower) half. You can also use your own matrix container, the error propagator will only call `::operator()(int,int)` on the container.

It is often possible to neglect correlations in a calculation. If the inputs are measurements, then they are likely independent (uncorrelated).
However, when a calculation is done with uncertain inputs, the result will also be correlated to the inputs. If you then do some calculation involving the results of a previous calculation,
then correlation *must* be tracked correctly.

The correlation matrix is a simple way to specify the correlation of inputs for a single function call, but it does not scale well. If two variables `x` and `y` are correlated and we calculate a variable `z` from these,
then the correlation between `x` and `y` can be specified in a 2 x 2 matrix. If we then did another calculation involving all three variables, we would need a 3 x 3 matrix to specify the correlation. We would also have
to be careful and ensure that the rows and columns of the matrix matches the argument number of each function we propagate error through.

libUncertainty provides provides a simple framework to help. To *automatically* tracking correlations between uncertain quantities, you can create a "correlation store", which is a container that stores the correlation between
pairs of uncertain variables. However, the plain `uncertain<...>` variables used in the previous examples will not work with a correlation store, we have to add an id.
```
correlation_store<double> store;

add_id<uncertain<double>> x(4, 0.1);
add_id<uncertain<double>> y(3, 0.2);

store.set(x,y,-1);

auto z = basic_error_propagator::propagate_error([](double a, double b) { return a + b; }, store, x, y);

z.nominal();      // 7
z.uncertainty();  // 0.1
store.get(z,x);   // -1
store.get(z,y);   // 1
```
In this example, `x` and `y` are uncertain variables with a unique ID. This ID is used by the store to track correlations. To do a calculation with correlated variables, we
add their correlation to the store, and then pass the store to the `propagate_error(...)` function. The function will automatically compute the correlation between the result and
each input, and add it to the store. If we start with uncorrelated inputs, then we would not need to add anything to the store.

Note that instances of `correlation_store<...>` are independent, they do not know about the correlations stored in other stores. If you want to use a common, shared correlation store for your entire application,
you can get a reference to a global store by calling `get_global_correlation_store()`.
```
auto& store = get_global_correlation_store();

add_id<uncertain<double>> x(4, 0.1);
add_id<uncertain<double>> y(3, 0.2);

store.set(x,y,-1);

auto z = basic_error_propagator::propagate_error([](double a, double b) { return a + b; }, store, x, y);

z.nominal();      // 7
z.uncertainty();  // 0.1
store.get(z,x);   // -1
store.get(z,y);   // 1
```

## Error Propagation Method

The library provides a simple error propagation method that is described in "An Introduction to Error Analysis" by John R. Taylor. It is a simple method that can be
done "by hand" (in fact, it is the method we teach students in our physics laboratory courses at Fort Hays State University), and it is equivalent to the first-order
theory for small deviations. See the [technical write up](./doc/writups/TechnicalDetails/TechnicalDetails.pdf) for details.
