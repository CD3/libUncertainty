# libUncertainty

A C++ library for working with uncertain quantities.

This is a _simple_ library, largely based on my Python library [pyErrorProp](https://github.com/CD3/pyErrorProp)
that provides some _basic_ utilities for working with uncertain quantities, including _basic_ error propagation.
It was written with [Boost.Units](https://www.boost.org/doc/libs/1_78_0/doc/html/boost_units.html) in mind, so
it supports uncertain quantities with units.

## Installation

Build and install with CMake:

    $ mkdir build
    $ cd build
    $ cmake ..
    $ cmake --build .
    $ cmake --build . --target install



## Usage

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
Calculations with uncertain variables is not directly supported. Instead, error propagation (through arbitrary functions)
is done with a wrapper function.
```
#include <libUncertainty/propagate.hpp>
...
uncertain<double> x(2, 0.2);
auto              y = error_propagators::basic::propagate_error([](double x) { return x * x; }, x);
y.nominal();     // 4
y.uncertainty(); // 0.84
```
This works with quantities too
```
uncertain<quantity<t::cm>> L(2 * i::cm, 0.1 * i::cm);
uncertain<quantity<t::cm>> W(4 * i::cm, 0.2 * i::cm);

auto A = error_propagators::basic::propagate_error([](quantity<t::cm> L, quantity<t::cm> W) { return L * W; }, L, W);

A.nominal();      // 8 cm
A.uncertainty();  // 0.5656754 cm
```

Currently, error propagation is very limited. Only functions with three or fewer arguments are supported, and correlated
inputs are NOT supported.
