To view these fractals, please visit [Mandelbrot-Qt](https://https://github.com/calum74/mandelbrot-qt) for an application which renders these fractals. This contains the internal data structures which may be of interest.

# Why a C++ template library

Although C++ templates are ugly and hard to read, they offer a lot of benefits.

A C++ template library allows the same code to work with many different configurations and data-types. If allows you to for example specify the power of the mandelbrot set `N`, the number of terms to compute in a series `T`, and the type of complex number you want to use for your calculation `Complex`. This will generate efficient code without sacrificing performance.

```c++
template<int N> struct mandelbrot_calculation
{
  ...
  template <typename Complex, unsigned long T>
  static std::array<Complex, T>
  delta_terms(const Complex &z, const std::array<Complex, T> &previous) {
     ...
  }
};
```

C++ templates allow different datatypes to be used at different resolutions, for example using `high_precision_real<3>` at low resolutions, and `high_precision_real<6>` at the next resolution.

Importantly, templates enable fine-tuning and experiments. You can experiment with a new number type (e.g. `high_exponent_real`) without rewriting a lot of code. You can experiment with adding more terms to a Taylor series just by tuning a template parameter.

# Classes

## high_exponent_real<Value, Exponent>

Extends a `double` with an additional `Exponent` field, for situations where a standard `double` precision number is not sufficient (for example deltas).

## high_precision_real<N>

Arbitrary-precision fixed-point arithmetic. `N` refers (quite awkwardly) to the number of 64-bit integers used to store the number. This is not optimized for ridiculously large values of `N`.

## mandelbrot_calculation<N>

Represents the equations for calculating a Mandelbrot set. `N` refers to the "power" of the Mandelbrot set, which is traditionally 2, but extends to arbitrary integers >1 (and indeed non-integers but that is out of scope).

$`z -> z^N + c`$

# About this repo

This repository contains algorithms for computing the Mandelbrot set. Although the basic Mandelbrot set is quite easy to compute, more advanced algorithms are necessary for good performance at very deep zooms, such as using high-precision arithmetic, perturbation theory and Taylor series approximations.

The algorithms have been organised as a template library, to make it easy to modify the algorithms to for example add different data types, or to have different data types for different levels of precision.

In addition, Mandelbrot orbits have been implemented as C++ iterators.

# Mandeldrop

The Mandeldrop fractal is a variation of the Mandelbrot set, where each point is iterated according to `z -> z^2 + (1/c)`. This essentially inverts the image such that the black is on the outside instead of the image. However, this is a conformal mapping, so zooming in, the image looks the same as a regular Mandelbrot set.

In order to rotate the image into a droplet, the image is transformed further by a multiplication of `-i` which merely rotates the image by π/2. So the final algorithm becomes `z -> z^2 + (-i/c)`. (Detail: the image is often rendered with the y-axis going downwards, hence the minus.)

In order for this to work with perturbation theory, a little bit of mathing gives that we need to use the delta `id/(c*(c+d))`, where `c` is the reference orbit and `d` is the original delta. Then we can simply use the standard Mandelbrot algorithms to compute the image.

I have seen some images of a Mandeldrop fractal before, so this idea is not original, but rotation through π/2 and use of perturbation theory might be new.

# High precision perturbation

Original work by Calum Grant.

After a certain point, say ay 10^-320, regular perturbation implementations are no longer sufficient because `double` numbers are no longer able to hold even the deltas. `long double` numbers are available on certain platforms, but not all.

To fix this, we can scale each delta with a factor M, which is a very small number like 2^-640. Then we rewrite our perturbation equations with d = M.D and e = M.E.

This gives

```
(1)     E_(n+1) = 2.z_n.E_n + M.E_n^2 + D
```

The Taylor series for E_n in terms of D around a central orbit is given by

```
(2)     E_n = A_n.D + B_n.D^2 + C_n.D^3 + O(D^4)
```

We'll assume that M and D are both very small, so we can ignore O(D^4).  Substituting (2) into (1),

```
(3)     E_(n+1) = 2.z_n.(A_n.D + B_n.D^2 + C_n.D^3 + O(D^4)) + M.(A_n.D + B_n.D^2 + C_n.D^3 + O(D^4))^2 + D

                = (2.z_n.A_n + 1).D + (2.z_n.B_n + M.A_n^2).D^2 + (2.z_n.C_n + 2.M.A_n.B_n).D^3 + O(D^4)
```

The Taylor series coefficients of E satisfy the iterative relation

```
(4)     A_(n+1) = 2.Z_n.A_n + 1
        B_(n+1) = 2.z_n.B_n + M.A_n^2
        C_(n+1) = 2.z_n.C_n + 2.M.A_n.B_n
        D_(n+1) = 2.z_n.D_n + 2.M.A_n.D_n + 2.M.b_n.C_n
        ...
```

We can use the adapted equations (1) and (4) in our Mandelbrot calculations to compute with higher precision.
