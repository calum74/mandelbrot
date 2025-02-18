/*
  Additional utilities and algorithms for dealing with complex numbers.
*/

#pragma once
#include "real_number.hpp"
#include <complex>

namespace fractals {

  template <typename T>
  concept Complex = requires(T v) {
    // requires(typename T::value_type) {};
    { v.real() } -> std::same_as<typename T::value_type>;
    v.imag();
  };
  
template <typename T> T real_part(const std::complex<T> &c) { return c.real(); }

template <typename T> T imag_part(const std::complex<T> &c) { return c.imag(); }

template <Complex C> C square(const C &c) {
  return {real_part(c) * real_part(c) - imag_part(c) * imag_part(c),
          typename C::value_type(2) * real_part(c) * imag_part(c)};
}

template <typename T>
std::complex<T> operator*(int k, const std::complex<T> &c) {
  return {T(k) * real_part(c), T(k) * imag_part(c)};
}

template <typename T>
std::complex<T> mul(const std::complex<T> &a, const std::complex<T> &b) {
  return {real_part(a) * real_part(b) - imag_part(a) * imag_part(b),
          imag_part(a) * real_part(b) + real_part(a) * imag_part(b)};
}

inline double to_double(double d) { return d; }
inline double to_double(float f) { return f; }

// template <Complex C> C step(const C &z, const C &c) { return square(z) + c; }

template <Complex C> auto norm(const C &c) {
  auto r = real_part(c);
  auto i = imag_part(c);
  return r * r + i * i;
}

template <int Order, Complex C, bool is_even = (Order % 2 == 0)>
struct pow_impl;

template <Complex C> struct pow_impl<2, C, true> {
  static C eval(const C &c) { return square(c); }
};

template <Complex C> struct pow_impl<1, C, false> {
  static C eval(const C &c) { return c; }
};

template <int Order, Complex C> struct pow_impl<Order, C, true> {
  static C eval(const C &c) {
    auto r = pow_impl<Order / 2, C>::eval(c);
    return square(r);
  }
};

template <Complex C> struct pow_impl<0, C, true> {
  static C eval(const C &c) { return C{1, 0}; }
};

template <Complex C> struct pow_impl<-1, C, false> {
  static C eval(const C &c) { return C{0, 0}; } // Not implemented
};

template <int Order, Complex C> struct pow_impl<Order, C, false> {
  static C eval(const C &c) { return mul(c, pow_impl<Order - 1, C>::eval(c)); }
};

template <int Order, Complex C> C pow(const C &c) {
  return pow_impl<Order, C>::eval(c);
}

template <int N, int M> struct choose_impl {
  static const int value =
      choose_impl<N - 1, M - 1>::value + choose_impl<N - 1, M>::value;
};

template <int M> struct choose_impl<0, M> {
  static const int value = 0;
};

template <> struct choose_impl<0, 1> {
  static const int value = 0;
};

template <int N> struct choose_impl<N, N> {
  static const int value = 1;
};

template <int N> struct choose_impl<N, 0> {
  static const int value = 1;
};

template <int N> struct choose_impl<N, 1> {
  static const int value = N;
};

template <int N, int M> constexpr int choose() {
  return choose_impl<N, M>::value;
}

inline bool isfinite(double d) { return std::isfinite(d); }

template <typename T> std::complex<T> normalize(const std::complex<T> &c) {
  return {normalize(c.real()), normalize(c.imag())};
}

template <typename T> struct normalized<std::complex<T>> {
  using type = std::complex<typename normalized<T>::type>;
};

// A complex number with the given precision
template <int Digits, int MinExp, int MaxExp>
using complex_number =
    std::complex<typename make_real<Digits, MinExp, MaxExp>::type>;


} // namespace fractals
