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

template <Complex C> bool escaped(const C &c) {
  return fractals::norm(c) >= typename C::value_type(4);
}

template <int N> struct mandelbrot_calculation;

namespace detail {
template <int N, int J> struct calculate_epsilon {
  template <Complex OrbitType, Complex DeltaType>
  static DeltaType eval(const OrbitType &z, const DeltaType &e) {
    return choose<N, J>() * DeltaType(pow<J>(z)) * pow<N - J>(e) +
           calculate_epsilon<N, J + 1>::eval(z, e);
  }
};

template <int N> struct calculate_epsilon<N, N> {
  template <Complex OrbitType, Complex DeltaType>
  static DeltaType eval(const OrbitType &z, const DeltaType &e) {
    return {0, 0};
  }
};

// We need to pick exactly one term from the equation
// (A∂ + B∂^2 + C∂^3 ...)
// and multiply it with the terms from (A∂ + B∂^2 + C∂^3
// ...)^(seq_remaining-1)
template <Complex TermType, unsigned long Terms>
void distribute_terms(int delta_pow, int seq_remaining, TermType f,
                      const std::array<TermType, Terms> &previous,
                      std::array<TermType, Terms> &next) {
  if (seq_remaining == 0) {
    next[delta_pow - 1] += f;
  } else {
    for (int p = 1; p + delta_pow + seq_remaining <= Terms + 1; p++) {
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

template <Complex TermType, unsigned long Terms, int N>
struct calculate_delta_terms {
  static std::array<TermType, Terms>
  calculate(TermType z, const std::array<TermType, Terms> &previous) {
    std::array<TermType, Terms> result;
    result[0] = 1;    // Represents the first ∂ in the equation for epsilon'
    TermType zJ(1, 0); // z^j

    for (int j = 0; j < N; j++, zJ = zJ * z) {
      int nCj = fac(N) / (fac(j) * fac(N - j));
      distribute_terms(0, N - j, zJ * TermType(nCj), previous, result);
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

// Specialisation of the previous case for N=4 (might be a bit faster?)
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
  template <Complex C> static C step(const C &z, const C &c) {
    return pow<N>(z) + c;
  }

  // When performing perturbations (for higher precision), here is the general
  // formula for evaluating the epsilon (dz)
  template <Complex OrbitType, Complex DeltaType>
  static DeltaType step_epsilon(const OrbitType &z, const DeltaType &e,
                                const DeltaType &d) {
    return d + detail::calculate_epsilon<N, 0>::eval(z, e);
  }

  template <Complex C> static C A(const C &z, const C &A_prev) {
    return choose<N, 1>() * pow<N - 1>(z) * A_prev + C{1, 0};
  }

  template <Complex C>
  static C B(const C &z, const C &A_prev, const C &B_prev) {
    return choose<N, 1>() * pow<N - 1>(z) * B_prev +
           choose<N, 2>() * pow<N - 2>(z) * pow<2>(A_prev);
  }

  template <Complex Cx>
  static Cx C(const Cx &z, const Cx &A_prev, const Cx &B_prev,
              const Cx &C_prev) {
    if constexpr (N > 2) {
      return choose<N, 1>() * pow<N - 1>(z) * C_prev +
             2 * choose<N, 2>() * pow<N - 2>(z) * A_prev * B_prev +
             choose<N, 3>() * pow<N - 3>(z) * pow<3>(A_prev);
    } else {
      return choose<N, 1>() * pow<N - 1>(z) * C_prev +
             2 * choose<N, 2>() * pow<N - 2>(z) * A_prev * B_prev;
    }
  }

  template <Complex C>
  static C D(const C &z, const C &A_prev, const C &B_prev, const C &C_prev,
             const C &D_prev) {
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
  //             = ∂ + sum(j=0..n-1)( choose(n,j).z^j.(A∂ + B∂^2 + C∂^3
  //             ...)^(n-j) )
  //
  // for the next iteration. We equate terms in ∂^n.
  template <Complex OrbitType, Complex TermType, unsigned long Terms>
  static std::array<TermType, Terms>
  delta_terms(const OrbitType &z, const std::array<TermType, Terms> &previous) {
    return detail::calculate_delta_terms<TermType, Terms, N>::calculate(
        TermType(z), previous);
  }
}; // mandelbrot_calculation

} // namespace mandelbrot