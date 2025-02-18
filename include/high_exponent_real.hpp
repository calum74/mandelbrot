#pragma once

#include <cmath>
#include <iostream>
#include <complex>

namespace fractals {

template <std::floating_point Double, std::integral Exponent, bool Normalize>
class high_exponent_real;

using high_exponent_double = high_exponent_real<double, int, false>;

using fast_high_exponent_double = high_exponent_real<double, int, false>;

/*
  An exponented real is a normal floating point number (like a double) with
  extra precision for the exponent. This is a good compromise when we need a bit
  more precision and still get good performance.

  Can in theory be combined with `high_precision_real` for an even higher
  number.
*/
template <std::floating_point Double = double, std::integral Exponent = int,
          bool Normalize = false>
class high_exponent_real {
public:
  template <bool FromNormalized>
  constexpr high_exponent_real(
      const high_exponent_real<Double, Exponent, FromNormalized> &src) {
    if constexpr (Normalize && !FromNormalized) {
      if (src.d == 0) {
        d = 0;
        e = 0;
      } else {
        d = std::frexp(src.d, &e);
        e += src.e;
      }
    } else {
      d = src.d;
      e = src.e;
    }
  }

  high_exponent_real(Double d0 = 0, Exponent e0 = 0) {
    if constexpr (Normalize) {
      if (d0 == 0) {
        d = 0;
        e = 0;
      } else {
        int e1;
        d = std::frexp(d0, &e1);
        e = e0 + e1;
      }
    } else {
      d = d0;
      e = e0;
    }
  }
  // explicit high_exponent_real(int i) : high_exponent_real(i, 0) {}

  Double d;
  Exponent e;

  Double to_double() const { return std::ldexp(d, e); }

  high_exponent_real &operator+=(high_exponent_real a) {
    return *this = (*this + a);
  }

  high_exponent_real &operator-=(high_exponent_real a) {
    return *this = (*this - a);
  }

  high_exponent_real &operator*=(high_exponent_real a) {
    return *this = (*this * a);
  }
};

template <typename D, typename E, bool N>
high_exponent_real<D, E, true> normalize(const high_exponent_real<D, E, N> x) {
  return x;
}

template <typename D, typename E, bool N>
high_exponent_real<D, E, false> operator-(high_exponent_real<D, E, N> a) {
  return high_exponent_real<D, E, false>{-a.d, a.e};
}

template <typename D, typename E, bool N1, bool N2>
high_exponent_real<D, E, false> operator-(high_exponent_real<D, E, N1> a,
                                          high_exponent_real<D, E, N2> b) {
  return a + -b;
}

template <typename D, typename E, bool N1, bool N2>
high_exponent_real<D, E, false> operator+(high_exponent_real<D, E, N1> a,
                                          high_exponent_real<D, E, N2> b) {
  if (a.d == 0)
    return b;
  if (b.d == 0)
    return a;

  // Need to figure out a sensible exponent
  if (a.e > b.e) {
    // Convert b to a's exponent
    return high_exponent_real<D, E, false>{a.d + b.d * (D)std::exp2(b.e - a.e),
                                           a.e};
  } else if (a.e < b.e) {
    // Convert a to b's exponent
    return high_exponent_real<D, E, false>{a.d * (D)std::exp2(a.e - b.e) + b.d,
                                           b.e};
  } else {
    return high_exponent_real<D, E, false>{a.d + b.d, a.e};
  }
}

template <typename D, typename E, bool N1, bool N2>
high_exponent_real<D, E, false> operator*(high_exponent_real<D, E, N1> a,
                                          high_exponent_real<D, E, N2> b) {
  return high_exponent_real<D, E, false>{a.d * b.d, a.e + b.e};
}

template <typename D, typename E, bool N1, bool N2>
high_exponent_real<D, E, false> operator/(high_exponent_real<D, E, N1> a,
                                          high_exponent_real<D, E, N2> b) {
  return high_exponent_real<D, E, false>{a.d / b.d, a.e - b.e};
}

template <typename D, typename E>
std::ostream &operator<<(std::ostream &os, high_exponent_real<D, E> a) {
  os << a.d;
  if (a.e)
    os << "*2^" << a.e;
  return os;
}

template <typename D, typename E>
std::istream &operator>>(std::istream &is, high_exponent_real<D, E> &a) {
  // Read up to the 'e'.
  bool negative;
  bool hexadecimal;

  // Read sign
  char p = is.peek();
  if (p == '-') {
    negative = true;
    is.get(p);
  } else {
    negative = false;
  }

  if (p == '0') {
    is.get(p);
    p = is.peek();
    if (p == 'x')
      hexadecimal = true;
    else
      hexadecimal = false;
  } else
    hexadecimal = false;

  // Keep reading until the '.'

  return is;
}

template <typename D, typename E>
int cmp(high_exponent_real<D, E, true> a, high_exponent_real<D, E, true> b) {
  // If the sign is different
  if (a.d < 0 && b.d > 0)
    return -1;
  if (a.d > 0 && b.d < 0)
    return 1;

  int a_bigger = a.d < 0 ? -1 : 1;

  // !! This does not work for zero !! FIXME
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

template <typename D, typename E, bool N1, bool N2>
bool operator==(high_exponent_real<D, E, N1> a,
                high_exponent_real<D, E, N2> b) {
  return cmp(normalize(a), normalize(b)) == 0;
}

template <typename D, typename E, bool N1, bool N2>
bool operator>=(high_exponent_real<D, E, N1> a,
                high_exponent_real<D, E, N2> b) {
  return cmp(normalize(a), normalize(b)) >= 0;
}

template <typename D, typename E, bool N1, bool N2>
bool operator>(high_exponent_real<D, E, N1> a, high_exponent_real<D, E, N2> b) {
  return cmp(normalize(a), normalize(b)) > 0;
}

template <typename D, typename E, bool N1, bool N2>
bool operator<(high_exponent_real<D, E, N1> a, high_exponent_real<D, E, N2> b) {
  return cmp(normalize(a), normalize(b)) < 0;
}

template <typename D, typename E, bool N1, bool N2>
bool operator<=(high_exponent_real<D, E, N1> a,
                high_exponent_real<D, E, N2> b) {
  return cmp(normalize(a), normalize(b)) <= 0;
}

template <typename D, typename E> bool isfinite(high_exponent_real<D, E> a) {
  return std::isfinite(a.d);
}

template <typename D, typename E, bool N>
D log(const high_exponent_real<D, E, N> &a) {
  return std::log(a.d) + a.e * std::log(2);
}

template <typename D, typename E, bool N>
high_exponent_real<D, E, false> operator<<(const high_exponent_real<D, E, N> &a,
                                           int shift) {
  return {a.d, a.e + shift};
}

template <typename D, typename E, bool N>
high_exponent_real<D, E, false> operator>>(const high_exponent_real<D, E, N> &a,
                                    int shift) {
  return a << -shift;
}


// Specialise this because std::complex does incompatible things
template <typename D, typename E, bool N1, bool N2>
std::complex<high_exponent_real<D, E, false>>
operator*(std::complex<high_exponent_real<D, E, N1>> a,
          std::complex<high_exponent_real<D, E, N2>> b) {
  return {a.real() * b.real() - a.imag() * b.imag(),
          a.real() * b.imag() + a.imag() * b.real()};
}

template <typename D, typename E, bool N1, bool N2>
std::complex<high_exponent_real<D, E, false>>
operator+(std::complex<high_exponent_real<D, E, N1>> a,
          std::complex<high_exponent_real<D, E, N2>> b) {
  return {a.real() + b.real(), a.imag() + b.imag()};
}

template <typename D, typename E, bool N>
std::complex<high_exponent_real<D, E, false>>
operator+(std::complex<high_exponent_real<D, E, N>> a,
          std::complex<high_exponent_real<D, E, N>> b) {
  return {a.real() + b.real(), a.imag() + b.imag()};
}

template <typename D, typename E, bool N>
std::complex<high_exponent_real<D, E, false>>
operator*(std::complex<high_exponent_real<D, E, N>> a,
          std::complex<high_exponent_real<D, E, N>> b) {
  return {a.real() * b.real() - a.imag() * b.imag(),
          a.real() * b.imag() + a.imag() * b.real()};
}

} // namespace fractals
