To view these fractals, please visit [Mandelbrot-Qt](https://https://github.com/calum74/mandelbrot-qt) for an application which renders these fractals. This repository contains the internal data structures which may be of interest.

# Mandelbrot C++ template library

## Overview

This library contains a variety of data structures and algorithms implemented using C++ templates.

Although C++ templates are harder to read, they give the additional flexibility needed to plug in new numerical types and allow for much better code reuse.

Unfortunately, this documentation is very superficial!

## List of classes

### complex_number.hpp

| Concept | Description |
| ----- | ----------- |
| `Complex` | Any complex number, for example `std::complex`. |

| Function | Description |
| ----- | ----------- |
| `pow<int N>()` | Raises a complex number to a fixed integer power in $lg N$ multiplications. |
| `choose<int N, int M>()` | Statically calculates ${N\choose M}$. |

| Type | Description |
| ----- | ----------- |
| `complex_number<int Bits, int MinExp, int MaxExp>` | A `std::complex` number type with the required precision and exponent. Very high precisions and exponents are supported. |

### fractal.hpp

Type | Description
-- | --
`pointwise_fractal` | A fractal.
`pointwise_calculation` | The logic to calculate a fractal.

Function | Description 
-- | --
`make_fractal()` | Creates a `pointwise_fractal`. 

### high_exponent_real.hpp

Type | Description
-- | --
`high_exponent_real<>` | A number that has fairly low precision, but supports very high exponents.
`high_exponent_double` | A `double` with a 32-bit exponent. Used for series terms which can get pretty large.

### high_precision_real.hpp

Type | Description
-- | --
`high_precision_real<int FractionalBits>` | Represents a real number to a given precision. Needed to calculate the reference orbit and to represent view coordinates.

### mandelbrot_calculation.hpp

| Type | Description |
| ----- | ----------- |
| `mandelbrot_calculation<int Power>` | Implements common formulae and perturbation expressions for a Mandelbrot set of the given power. |

| Function | Description |
| ----- | ----------- |
| `escaped()` | Tests for escape. |
| `mandelbrot_calculation::step()` | Performs a single iteration. |
| `mandelbrot_calculation::step_epsilon()` | Computes $\epsilon_{i+1}$ ($\Delta z_{i+1}$). |
| `mandelbrot_calculation::delta_terms()` | Computes the Taylor series terms to arbitrary length for any power of Mandelbrot set. |

### orbit.hpp

| Concept | Description |
| ----- | ----------- |
| `Calculation` | A type providing Mandelbrot term calculations. | 
| `IteratedOrbit` | An orbit supporting `++`. |
| `RandomAccessOrbit` | An orbit supporting `[]`. |

| Type | Description |
| ----- | ----------- |
| `basic_orbit<>` | An `IteratedOrbit`, calculated using the standard formula, typically used for the high precision reference orbit. | 
| `stored_orbit<>` | A `RandomAccessOrbit`, stores all of the values from a `basic_orbit<>` in an array, typically using low precision numbers. |
| `taylor_series_orbit<>` | An `IteratedOrbit` that also calculates the Taylor series terms to any length. |
| `perturbation_orbit<>` | An `IteratedOrbit` relative to a `stored_orbit<>`. Is able to reset the orbit which is why a `stored_orbit<>` is needed as the reference orbit. |
| `stored_taylor_series_orbit<>` | A `RandomAccessOrbit` which stores the Taylor series terms in an array, and is able to skip forward to create a `perturbation_orbit<>`. |

### percentile.hpp

Function | Description
-- | --
`top_percentile()` | Returns an element from an unsorted array at the given percentile, more efficiently than by sorting it.

### plane.hpp

Type | Description
-- | --
`plane<>` | The mapping between pixel coordinates and complex coordinates.

### real_number.hpp

Type | Description
-- | --
`real_number<int Bits, int MinExp, int MaxExp>` | A numerical datatype that is able to represent the specified level of precision.

### view_coords.hpp

Type | Description
-- | --
`view_coords` | The size and position of a view.

### view_parameters.hpp

Type | Description
-- | --
`view_parameters` | All parameters required to specify a view.
