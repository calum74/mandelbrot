#pragma once
#include "convert.hpp"
#include "high_precision_real.hpp"
#include <complex>

#include <cmath>
#include <iostream>

namespace fractals {

struct skip_normalization {};
/*
  An exponented real is a normal floating point number (like a double) with
  extra precision for the exponent. This is a good compromise when we need a bit
  more precision and still get good performance.

  Can in theory be combined with `high_precision_real` for an even higher
  number.
*/
template <typename Double = double, typename Exponent = int>
class high_exponent_real {
public:
  high_exponent_real(Double d0, Exponent e0, skip_normalization)
      : d(d0), e(e0) {}

  high_exponent_real(Double d0 = 0, Exponent e0 = 0) {
    if (d0 == 0) {
      d = 0;
      e = 0;
    } else {
      int e1;
      d = std::frexp(d0, &e1);
      e = e0 + e1;
    }
  }
  high_exponent_real(int i) : high_exponent_real(i, 0) {}

  Double d;
  Exponent e;

  double to_double() const { return std::ldexp(d, e); }

  high_exponent_real &operator+=(high_exponent_real a) {
    return *this = (*this + a);
  }
};

template <typename D, typename E>
high_exponent_real<D, E> operator-(high_exponent_real<D, E> a) {
  return {-a.d, a.e};
}

template <typename D, typename E>
high_exponent_real<D, E> operator-(high_exponent_real<D, E> a,
                                   high_exponent_real<D, E> b) {
  return a + -b;
}

template <typename D, typename E>
high_exponent_real<D, E> operator+(high_exponent_real<D, E> a,
                                   high_exponent_real<D, E> b) {
  // Need to figure out a sensible exponent
  if (a.e > b.e) {
    // Convert b to a's exponent
    return {a.d + b.d * std::exp2(b.e - a.e), a.e};
  } else if (a.e < b.e) {
    // Convert a to b's exponent
    return {a.d * std::exp2(a.e - b.e) + b.d, b.e};
  } else {
    return {a.d + b.d, a.e};
  }
}

template <typename D, typename E>
high_exponent_real<D, E> operator*(high_exponent_real<D, E> a,
                                   high_exponent_real<D, E> b) {
  return {a.d * b.d, a.e + b.e};
}

template <typename D, typename E>
high_exponent_real<D, E> operator/(high_exponent_real<D, E> a,
                                   high_exponent_real<D, E> b) {
  return high_exponent_real<D, E>{a.d / b.d, a.e - b.e};
}

template <typename D, typename E>
std::ostream &operator<<(std::ostream &os, high_exponent_real<D, E> a) {
  os << a.d;
  if (a.e)
    os << "*2^" << a.e;
  return os;
}

template <typename D, typename E>
int cmp(high_exponent_real<D, E> a, high_exponent_real<D, E> b) {

  // If the sign is different
  if (a.d < 0 && b.d > 0)
    return -1;
  if (a.d > 0 && b.d < 0)
    return 1;

  int a_bigger = a.d < 0 ? -1 : 1;

  // They are the same sign
  // Compare exponents
  if (a.e > b.e)
    return a_bigger;
  if (a.e < b.e)
    return -a_bigger;

  // The exponent is the same
  if (a.d < b.d)
    return -1;
  if (a.d > b.d)
    return 1;
  return 0;
}

template <typename D, typename E>
bool operator==(high_exponent_real<D, E> a, high_exponent_real<D, E> b) {
  return a.d == b.d && a.e == b.e;
}

template <typename D, typename E>
bool operator>=(high_exponent_real<D, E> a, high_exponent_real<D, E> b) {
  return cmp(a, b) >= 0;
}

template <typename D, typename E>
bool operator>(high_exponent_real<D, E> a, high_exponent_real<D, E> b) {
  return cmp(a, b) > 0;
}

template <typename D, typename E>
bool operator<(high_exponent_real<D, E> a, high_exponent_real<D, E> b) {
  return cmp(a, b) < 0;
}

template <typename D, typename E>
bool operator<=(high_exponent_real<D, E> a, high_exponent_real<D, E> b) {
  return cmp(a, b) <= 0;
}

template <typename D, typename E>
struct convert_to<D, high_exponent_real<D, E>> {
  static D get(const high_exponent_real<D, E> &x) { return x.to_double(); }
};

template <int N, typename D, typename E>
struct convert_to<high_precision_real<N>, high_exponent_real<D, E>> {
  static high_precision_real<N> get(const high_exponent_real<D, E> &x) {
    return high_precision_real<N>(x.d) << x.e;
  }
};

template <int N, typename D, typename E>
struct convert_to<high_exponent_real<D, E>, high_precision_real<N>> {
  static high_exponent_real<D, E> get(const high_precision_real<N> &x) {
    if (x.fraction[0] == 0) {
      int e = count_fractional_zeros(x);
      return {(x >> e).to_double(), e};
    } else {
      return x.to_double();
    }
  }
};

template <typename D, typename E> bool isfinite(high_exponent_real<D, E> a) {
  return std::isfinite(a.d);
}

template <typename D, typename E> D log(const high_exponent_real<D, E> &a) {
  return std::log(a.d) + a.e * std::log(2); // Test this !!
}

template <typename D, typename E>
high_exponent_real<D, E> operator<<(const high_exponent_real<D, E> &a,
                                    int shift) {
  return {a.d, a.e + shift};
}

template <typename D, typename E>
high_exponent_real<D, E> operator>>(const high_exponent_real<D, E> &a,
                                    int shift) {
  return a << -shift;
}

} // namespace fractals

// Override this because std::complex does incompatible things
template <typename D, typename E>
std::complex<fractals::high_exponent_real<D, E>>
operator*(std::complex<fractals::high_exponent_real<D, E>> a,
          std::complex<fractals::high_exponent_real<D, E>> b) {
  return {a.real() * b.real() - a.imag() * b.imag(),
          a.real() * b.imag() + a.imag() * b.real()};
}
