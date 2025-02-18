#pragma once
#include "complex_number.hpp"

namespace fractals {
template <typename To, typename From> struct convert_to;

template <typename T> struct convert_to<T, T> {
  static T get(const T &t) { return t; }
};

template <> struct convert_to<double, int> {
  static double get(int i) { return i; }
};

template <> struct convert_to<double, float> {
  static double get(float f) { return f; }
};

template <> struct convert_to<float, double> {
  static float get(double f) { return f; }
};

template <typename T1, typename T2> T1 convert(const T2 &x) {
  return convert_to<T1, T2>::get(x);
}

template <int N, int M>
struct convert_to<high_precision_real<N>, high_precision_real<M>> {
  static high_precision_real<N> get(const high_precision_real<M> &x) {
    return x;
  }
};

template <int N>
struct convert_to<high_precision_real<N>, high_precision_real<N>> {
  static high_precision_real<N> get(const high_precision_real<N> &x) {
    return x;
  }
};

template <int M> struct convert_to<double, high_precision_real<M>> {
  static double get(const high_precision_real<M> &x) { return x.to_double(); }
};

template <int M> struct convert_to<float, high_precision_real<M>> {
  static float get(const high_precision_real<M> &x) { return x.to_double(); }
};

template <int N> struct convert_to<high_precision_real<N>, int> {
  static high_precision_real<N> get(int x) { return {x}; }
};

template <int M> struct convert_to<high_precision_real<M>, double> {
  static high_precision_real<M> get(double x) { return {x}; }
};

template <int M> struct convert_to<high_precision_real<M>, float> {
  static high_precision_real<M> get(float x) { return {x}; }
};

template <typename T1, typename T2>
  requires(!std::same_as<T1, T2>)
struct convert_to<std::complex<T1>, std::complex<T2>> {
  static std::complex<T1> get(const std::complex<T2> &c) {
    return {convert<T1>(c.real()), convert<T1>(c.imag())};
  }
};

template <std::floating_point D, typename E, bool N1, bool N2>
  requires(N1 != N2)
struct convert_to<high_exponent_real<D, E, N1>, high_exponent_real<D, E, N2>> {
  static high_exponent_real<D, E, N1> get(high_exponent_real<D, E, N2> x) {
    return x;
  }
};

template <std::floating_point D1, std::floating_point D2, typename E, bool N>
struct convert_to<D1, high_exponent_real<D2, E, N>> {
  static D1 get(const high_exponent_real<D2, E, N> &x) { return x.to_double(); }
};

template <std::floating_point D1, std::floating_point D2, typename E, bool N>
struct convert_to<high_exponent_real<D1, E, N>, D2> {
  static high_exponent_real<D1, E, N> get(D2 x) {
    return high_exponent_real<D1, E, N>{(D1)x};
  }
};

template <typename D, typename E, bool N>
struct convert_to<high_exponent_real<D, E, N>, int> {
  static high_exponent_real<D, E, N> get(int x) { return {x}; }
};

template <int N, typename D, typename E, bool Norm>
struct convert_to<high_precision_real<N>, high_exponent_real<D, E, Norm>> {
  static high_precision_real<N> get(const high_exponent_real<D, E, Norm> &x) {
    return high_precision_real<N>(x.d) << x.e;
  }
};

template <int N, std::floating_point D, typename E, bool Norm>
struct convert_to<high_exponent_real<D, E, Norm>, high_precision_real<N>> {
  static high_exponent_real<D, E, Norm> get(const high_precision_real<N> &x) {
    if (x.fraction[0] == 0) {
      int e = count_fractional_zeros(x);
      return high_exponent_real<D, E, Norm>{(D)(x << e).to_double(), -e};
    } else {
      return high_exponent_real<D, E, Norm>{x.to_double()};
    }
  }
};


} // namespace fractals