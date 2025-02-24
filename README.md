To view these fractals, please visit [Mandelbrot-Qt](https://https://github.com/calum74/mandelbrot-qt) for an application which renders these fractals. This repository contains the internal data structures which may be of interest.

# Mandelbrot C++ template library

## Overview

This library contains a variety of data structures and algorithms implemented using C++ templates.

Although C++ templates are harder to read, they give the additional flexibility needed to plug in new numerical types and allow for much better code reuse.

Unfortunately, this documentation is very superficial!

To see this library in action, see [tutorial.cpp](test/tutorial.cpp). You can also implement your own fractals, see the sample [circle.cpp](test/circle.cpp).

## List of concepts

| Concept | Header | Description |
| ----- | -- | ----------- |
| `Calculation` | orbit.hpp | A type providing Mandelbrot term calculations. | 
| `Complex` | complex_number.hpp | Any complex number, for example `std::complex`. |
| `IteratedOrbit` | orbit.hpp | An orbit supporting `++`. |
| `RandomAccessOrbit` | orbit.hpp | An orbit supporting `[]`. |

## List of types

| Type | Header | Description |
| ----- | -- | ----------- |
| `basic_orbit<>` | orbit.hpp | An `IteratedOrbit`, calculated using the standard formula, typically used for the high precision reference orbit. | 
| `complex_number<int Bits, int MinExp, int MaxExp>` | complex_number.hpp | A `std::complex` number type with the required precision and exponent. Very high precisions and exponents are supported. |
`fractal` | fractal.hpp | A fractal.
`fractal_calculation` | fractal_calculation.hpp | The logic to calculate the points of a fractal.
`high_exponent_double` | high_exponent_real.hpp | A `double` with a 32-bit exponent. Used for series terms which can get pretty large.
`high_exponent_real<>` | high_exponent_real.hpp | A number that has fairly low precision, but supports very high exponents.
`high_precision_real<int FractionalBits>` | high_precision_real.hpp | A real number to a given precision. Needed to calculate the reference orbit and to store view coordinates.
| `mandelbrot_calculation<int Power>` | mandelbrot_calculation.hpp | Implements common formulae and perturbation expressions for a Mandelbrot set of the given power. |
| `perturbation_orbit<>` | orbit.hpp | An `IteratedOrbit` relative to a `stored_orbit<>`. Is able to reset the orbit which is why a `stored_orbit<>` is needed as the reference orbit. |
`plane<>` | plane.hpp | The mapping between pixel coordinates and complex coordinates.
`real_number<int Bits, int MinExp, int MaxExp>` | real_number.hpp | A numerical datatype that is able to represent the specified level of precision.
| `stored_orbit<>` | orbit.hpp | A `RandomAccessOrbit`, stores all of the values from a `basic_orbit<>` in an array, typically using low precision numbers. |
| `taylor_series_orbit<>` | orbit.hpp | An `IteratedOrbit` that also calculates the Taylor series terms to any length. |
| `stored_taylor_series_orbit<>` | orbit.hpp | A `RandomAccessOrbit` which stores the Taylor series terms in an array, and is able to skip forward to create a `perturbation_orbit<>`. |
`view_coords` | view_coords.hpp | The size and position of a view.
`view_parameters` | view_parameters.hpp | All parameters required to specify a view.

## List of functions

| Function | Header | Description |
| -----| -- | ----------- |
| `choose<int N, int M>()` | complex_number.hpp | Statically calculates ${N\choose M}$. |
| `escaped()` | complex_number.hpp | Tests for escape. |
`make_fractal<...>()` | make_fractal.hpp | Creates a `fractal`. 
`number_cast<T>()` | number_cast.hpp | Converts a number from one type to another.
| `pow<int N>()` | complex_number.hpp | Raises a complex number to a fixed integer power in $lg N$ multiplications. |
`top_percentile()` | percentile.hpp | Returns an element from an unsorted array at the given percentile in $O(n)$ time.

