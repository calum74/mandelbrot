/*
  Additional utilities and algorithms for dealing with complex numbers.
*/

#pragma once
#include <complex>

namespace mandelbrot {

template <typename T> T real_part(const std::complex<T> &c) { return c.real(); }

template <typename T> T imag_part(const std::complex<T> &c) { return c.imag(); }

template <typename Complex> Complex square(const Complex &c) {
  return {real_part(c) * real_part(c) - imag_part(c) * imag_part(c),
          2 * real_part(c) * imag_part(c)};
}

template <typename T>
std::complex<T> operator*(int k, const std::complex<T> &c) {
  return {k * real_part(c), k * imag_part(c)};
}

template <typename T>
std::complex<T> mul(const std::complex<T> &a, const std::complex<T> &b) {
  return {real_part(a) * real_part(b) - imag_part(a) * imag_part(b),
          imag_part(a) * real_part(b) + real_part(a) * imag_part(b)};
}

inline double to_double(double d) { return d; }
inline double to_double(float f) { return f; }

template <typename T1> struct convert_to;

template <> struct convert_to<double> {
  template <typename T> static double get(const T &v) { return to_double(v); }
};

template <typename T1, typename T2> T1 convert_complex(const T2 &c) {
  return {convert_to<typename T1::value_type>::get(real_part(c)),
          convert_to<typename T1::value_type>::get(imag_part(c))};
}

template <typename C> C step(const C &z, const C &c) { return square(z) + c; }

template <typename C> auto norm(const C &c) {
  auto r = real_part(c);
  auto i = imag_part(c);
  return r * r + i * i;
}

template <typename C> bool escaped(const C &c) { return norm(c) >= 4; }

template <int Order, typename C, bool is_even = (Order % 2 == 0)>
struct pow_impl;

template <typename C> struct pow_impl<2, C, true> {
  static C eval(const C &c) { return square(c); }
};

template <typename C> struct pow_impl<1, C, false> {
  static C eval(const C &c) { return c; }
};

template <int Order, typename C> struct pow_impl<Order, C, true> {
  static C eval(const C &c) {
    auto r = pow_impl<Order / 2, C>::eval(c);
    return square(r);
  }
};

template <typename C> struct pow_impl<0, C, true> {
  static C eval(const C &c) { return C{1, 0}; }
};

template <typename C> struct pow_impl<-1, C, false> {
  static C eval(const C &c) { return C{0, 0}; } // Not implemented
};

template <int Order, typename C> struct pow_impl<Order, C, false> {
  static C eval(const C &c) { return mul(c, pow_impl<Order - 1, C>::eval(c)); }
};

template <int Order, typename C> C pow(const C &c) {
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

} // namespace mandelbrot
