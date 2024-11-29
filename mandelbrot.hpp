#pragma once

#include <cstddef>
#include <tuple>

namespace mandelbrot {

template <typename T> struct complex {
  T re={}, im={};
};

template <typename T> complex<T> square(complex<T> c) {
  return {c.re * c.re - c.im * c.im, 2 * c.re * c.im};
}

template <typename T> struct high_precision {
  T mantissa, exponent;
};

template <typename T>
high_precision<T> operator+(high_precision<T> a, high_precision<T> b) {
  auto e_delta = a.exponent - b.exponent;

  return e_delta < 0 ? high_precision<T>{(a.mantissa << e_delta) + b.mantissa, b.exponent}
                     : high_precision<T>{a.mantissa + (b.mantissa << e_delta), a.exponent};
}

template <typename T>
high_precision<T> operator*(high_precision<T> a, high_precision<T> b) {
  // What about overflows??
  return {a.mantissa * b.mantissa, a.exponent + a.exponent};
}

template <typename T> complex<T> operator*(complex<T> c1, complex<T> c2) {
  return {c1.re * c2.re - c1.im * c2.im, c1.im * c2.re + c1.re * c2.im};
}

template <typename T> complex<T> operator+(complex<T> c1, complex<T> c2) {
  return {c1.re + c2.re, c1.im + c2.im};
}

template<typename T> bool operator==(complex<T> c1, complex<T> c2)
{
    return std::tie(c1.re, c1.im) == std::tie(c2.re, c2.im);
}

template <typename C> C step(C z, C c) { return square(z) + c; }

template <typename T> bool escaped(complex<T> c);

template <typename C> int naive_count(C c, int max_iterations) {
  int i = 0;
  C z = c;
  while (++i < max_iterations && !escaped(z))
    z = step(z, c);
  return i;
}

} // namespace mandelbrot
