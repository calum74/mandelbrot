#pragma once
#include <complex>

namespace mandelbrot {

template <typename T> T real_part(const std::complex<T> &c) { return c.real(); }

template <typename T> T imag_part(const std::complex<T> &c) { return c.imag(); }

template <typename Complex> Complex square(const Complex &c) {
  return {real_part(c) * real_part(c) - imag_part(c) * imag_part(c),
          2 * real_part(c) * imag_part(c)};
}

template <typename T> struct complex_traits;

template <typename T> struct complex_traits<std::complex<T>> {
  using value_type = T;
};

template <typename T>
std::complex<T> operator*(int k, const std::complex<T> &c) {
  return {k * real_part(c), k * imag_part(c)};
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

} // namespace mandelbrot
