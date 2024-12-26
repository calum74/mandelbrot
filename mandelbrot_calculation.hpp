#pragma once

#include "complex.hpp"

namespace mandelbrot {

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

  template <int J> struct calculate_epsilon {
    template <typename Complex>
    static Complex eval(const Complex &z, const Complex &e) {
      return choose<N, J>() * pow<J>(z) * pow<N - J>(e) +
             calculate_epsilon<J + 1>::eval(z, e);
    }
  };

  template <> struct calculate_epsilon<N> {
    template <typename Complex>
    static Complex eval(const Complex &z, const Complex &e) {
      return {0, 0};
    }
  };

  // When performing perturbations (for higher precision), here is the general
  // formula for evaluating the epsilon (dz)
  template <typename Complex>
  static Complex step_epsilon(const Complex &z, const Complex &e,
                              const Complex &d) {
    return d + calculate_epsilon<0>::eval(z, e);
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
           choose<N, 2>() * pow<N - 2>(z) * Complex{2} * A_prev * B_prev +
           (N > 2 ? choose<N, 3>() * pow<N - 3>(z) * pow<3>(A_prev) : 0);
  }

  // TODO: Can generalize this further to calculate the Nth Taylor series term
};

} // namespace mandelbrot