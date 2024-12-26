/*
  Calculations involved in the Mandelbrot set.

  This is a "strategy pattern" that we'll plug into the various orbital
  algorithms in order to configure the algorithm, for example to vary the power
  of the Mandelbrot set.

  As well as the basic algorithm `z -> z^n + c`, we also include the calculation
  of deltas and Taylor series coefficeints for Mandelbrot sets of arbitrary
  powers.

  Only the first 3 Taylor series coefficients are calculated, in principle we
  could generalise this to an arbitrary number of coefficients although this has
  not been implemented.

  I have not seen these formulae derived anywhere but they are an obvious
  generalisation of the original work by superfractalthing by K.I. Martin.
  https://en.wikipedia.org/wiki/Plotting_algorithms_for_the_Mandelbrot_set
*/

#pragma once

#include "complex.hpp"

namespace mandelbrot {

namespace detail {
template <int N, int J> struct calculate_epsilon {
  template <typename Complex>
  static Complex eval(const Complex &z, const Complex &e) {
    return choose<N, J>() * pow<J>(z) * pow<N - J>(e) +
           calculate_epsilon<N, J + 1>::eval(z, e);
  }
};

template <int N> struct calculate_epsilon<N, N> {
  template <typename Complex>
  static Complex eval(const Complex &z, const Complex &e) {
    return {0, 0};
  }
};
} // namespace detail

// The complex arithmetic required to calculate a Mandelbrot set
// We consider the generalized Mandelbrot set, including higher orders
// such as the cubic, but exclude non-integer orders.
// TODO: Specialise this for N=2 for speed
template <int N> struct mandelbrot_calculation {

  static constexpr int order = N;

  // The general form of the calculation z -> z^N + c
  template <typename Complex>
  static Complex step(const Complex &z, const Complex &c) {
    return pow<N>(z) + c;
  }

  // When performing perturbations (for higher precision), here is the general
  // formula for evaluating the epsilon (dz)
  template <typename Complex>
  static Complex step_epsilon(const Complex &z, const Complex &e,
                              const Complex &d) {
    return d + detail::calculate_epsilon<N, 0>::eval(z, e);
  }

  template <typename Complex>
  static Complex A(const Complex &z, const Complex &A_prev) {
    return choose<N, 1>() * pow<N - 1>(z) * A_prev + Complex{1, 0};
  }

  template <typename Complex>
  static Complex B(const Complex &z, const Complex &A_prev,
                   const Complex &B_prev) {
    return choose<N, 1>() * pow<N - 1>(z) * B_prev +
           choose<N, 2>() * pow<N - 2>(z) * pow<2>(A_prev);
  }

  template <typename Complex>
  static Complex C(const Complex &z, const Complex &A_prev,
                   const Complex &B_prev, const Complex &C_prev) {
    return choose<N, 1>() * pow<N - 1>(z) * C_prev +
           2 * choose<N, 2>() * pow<N - 2>(z) * A_prev * B_prev +
           (N > 2 ? choose<N, 3>() * pow<N - 3>(z) * pow<3>(A_prev) : 0);
  }

  // TODO: Can generalize this further to calculate the Nth Taylor series term
};

} // namespace mandelbrot