/*
  Calculations involved in the Mandelbrot set.

  This is a "strategy pattern" that we'll plug into the various orbital
  algorithms in order to configure the algorithm, for example to vary the power
  of the Mandelbrot set.

  As well as the basic algorithm `z -> z^n + c`, we also include the calculation
  of deltas and Taylor series coefficeints for Mandelbrot sets of arbitrary
  powers.

  I have not seen these formulae derived anywhere but they are an obvious
  generalisation of the original work by superfractalthing by K.I. Martin.
  https://en.wikipedia.org/wiki/Plotting_algorithms_for_the_Mandelbrot_set
*/

#pragma once

#include "complex_number.hpp"
#include <array>

namespace mandelbrot {
using namespace fractals;

template <typename C> bool escaped(const C &c) {
  return fractals::norm(c) >= typename C::value_type(4);
}

template <int N> struct mandelbrot_calculation;

namespace detail {
template <int N, int J> struct calculate_epsilon {
  template <typename Complex, typename Delta>
  static Delta eval(const Complex &z, const Delta &e) {
    return choose<N, J>() * Delta(pow<J>(z)) * pow<N - J>(e) +
           calculate_epsilon<N, J + 1>::eval(z, e);
  }
};

template <int N> struct calculate_epsilon<N, N> {
  template <typename Complex, typename Delta>
  static Delta eval(const Complex &z, const Delta &e) {
    return {0, 0};
  }
};

// We need to pick exactly one term from the equation
// (A∂ + B∂^2 + C∂^3 ...)
// and multiply it with the terms from (A∂ + B∂^2 + C∂^3
// ...)^(seq_remaining-1)
template <typename Complex, unsigned long T>
void distribute_terms(int delta_pow, int seq_remaining, Complex f,
                      const std::array<Complex, T> &previous,
                      std::array<Complex, T> &next) {
  if (seq_remaining == 0) {
    next[delta_pow - 1] += f;
  } else {
    for (int p = 1; p + delta_pow + seq_remaining <= T + 1; p++) {
      // Distribute one term
      distribute_terms(delta_pow + p, seq_remaining - 1, f * previous[p - 1],
                       previous, next);
    }
  }
}

inline int fac(int n) {
  int m = 1;
  while (n > 1)
    m *= n--;
  return m;
}

template <typename Complex, unsigned long T, int N>
struct calculate_delta_terms {
  static std::array<Complex, T>
  calculate(Complex z, const std::array<Complex, T> &previous) {
    std::array<Complex, T> result;
    result[0] = 1;    // Represents the first ∂ in the equation for epsilon'
    Complex zJ(1, 0); // z^j

    for (int j = 0; j < N; j++, zJ = zJ * z) {
      int nCj = fac(N) / (fac(j) * fac(N - j));
      distribute_terms(0, N - j, zJ * Complex(nCj), previous, result);
    }

    return result;
  }
};

// Specialisation of the previous case for N=3 (might be a bit faster?)
template <typename Complex, int N> struct calculate_delta_terms<Complex, 3, N> {
  static std::array<Complex, 3>
  calculate(Complex z, const std::array<Complex, 3> &previous) {
    return {
        mandelbrot_calculation<N>::A(z, previous[0]),
        mandelbrot_calculation<N>::B(z, previous[0], previous[1]),
        mandelbrot_calculation<N>::C(z, previous[0], previous[1], previous[2])};
  }
};

// Specialisation of the previous case for N=3 (might be a bit faster?)
template <typename Complex, int N> struct calculate_delta_terms<Complex, 4, N> {
  static std::array<Complex, 4>
  calculate(Complex z, const std::array<Complex, 4> &previous) {
    return {
        mandelbrot_calculation<N>::A(z, previous[0]),
        mandelbrot_calculation<N>::B(z, previous[0], previous[1]),
        mandelbrot_calculation<N>::C(z, previous[0], previous[1], previous[2]),
        mandelbrot_calculation<N>::D(z, previous[0], previous[1], previous[2],
                                     previous[3])};
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
  template <typename Complex, typename Delta>
  static Delta step_epsilon(const Complex &z, const Delta &e, const Delta &d) {
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
    if constexpr (N > 2) {
      return choose<N, 1>() * pow<N - 1>(z) * C_prev +
             2 * choose<N, 2>() * pow<N - 2>(z) * A_prev * B_prev +
             choose<N, 3>() * pow<N - 3>(z) * pow<3>(A_prev);
    } else {
      return choose<N, 1>() * pow<N - 1>(z) * C_prev +
             2 * choose<N, 2>() * pow<N - 2>(z) * A_prev * B_prev;
    }
  }

  template <typename Complex>
  static Complex D(const Complex &z, const Complex &A_prev,
                   const Complex &B_prev, const Complex &C_prev,
                   const Complex &D_prev) {
    if constexpr (N > 3) {
      return choose<N, 1>() * pow<N - 1>(z) * D_prev +
             choose<N, 2>() * pow<N - 2>(z) *
                 (pow<2>(B_prev) + 2 * A_prev * C_prev) +
             choose<N, 4>() * pow<N - 4>(z) * pow<4>(A_prev);
    } else {
      return choose<N, 1>() * pow<N - 1>(z) * D_prev +
             choose<N, 2>() * pow<N - 2>(z) *
                 (pow<2>(B_prev) + 2 * A_prev * C_prev);
    }
  }

  // For a Taylor series expansion of the form
  //
  //    epsilon = A∂ + B∂^2 + C∂^3 ...
  //
  // compute the terms of
  //
  //    epsilon' = A'∂ + B'∂^2 + C'∂^3 ...
  //             = ∂ + sum(j=0..n-1)( C(n,j).z^j.(A∂ + B∂^2 + C∂^3 ...)^(n-j) )
  //
  // for the next iteration. We equate terms in ∂^n.
  template <typename Complex, unsigned long T>
  static std::array<Complex, T>
  delta_terms(const Complex &z, const std::array<Complex, T> &previous) {
    return detail::calculate_delta_terms<Complex, T, N>::calculate(z, previous);
  }
}; // mandelbrot_calculation

} // namespace mandelbrot