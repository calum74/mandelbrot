#pragma once

#include <cmath>
#include <complex>
#include <iostream>

namespace fractals {

/*
  An exponented real is a normal floating point number (like a double) with
  extra precision for the exponent. This is a good compromise when we need a bit
  more precision and still get good performance.
*/
template <std::floating_point Double = double, std::integral Exponent = int,
          bool Normalized = false>
class high_exponent_real {
public:
  template <bool SrcNormalized>
  constexpr high_exponent_real(
      const high_exponent_real<Double, Exponent, SrcNormalized> &src) {
    if constexpr (!SrcNormalized) {
      assign(src.mantissa(), src.exponent());
    } else {
      d = src.mantissa();
      e = src.exponent();
    }
  }

  template <bool SrcNormalized>
  high_exponent_real &
  operator=(const high_exponent_real<Double, Exponent, SrcNormalized> &src) {
    if constexpr (!SrcNormalized) {
      assign(src.mantissa(), src.exponent());
    } else {
      d = src.d;
      e = src.e;
    }
    return *this;
  }

  high_exponent_real<Double, Exponent, true> normalize() const { return *this; }

  void assign(Double d0, Exponent e0) {
    if constexpr (Normalized) {
      if (d0 == 0) { // !! Special case needed??
        d = 0;
        e = 0;
      } else {
        d = std::frexp(d0, &e);
        e += e0;
      }
    } else {
      d = d0;
      e = e0;
    }
  }

  high_exponent_real(Double d0 = 0, Exponent e0 = 0) { assign(d0, e0); }

  Double to_double() const { return std::ldexp(d, e); }
  Double mantissa() const { return d; }
  Exponent exponent() const { return e; }

  high_exponent_real &operator+=(high_exponent_real a) {
    return *this = (*this + a);
  }

  high_exponent_real &operator-=(high_exponent_real a) {
    return *this = (*this - a);
  }

  high_exponent_real &operator*=(high_exponent_real a) {
    return *this = (*this * a);
  }

private:
  Double d;
  Exponent e;
};

template <typename D, typename E, bool N>
high_exponent_real<D, E, true> normalize(const high_exponent_real<D, E, N> x) {
  return x;
}

template <typename D, typename E, bool N>
high_exponent_real<D, E, false> operator-(high_exponent_real<D, E, N> a) {
  return {-a.mantissa(), a.exponent()};
}

template <typename D, typename E, bool N1, bool N2>
high_exponent_real<D, E, false> operator-(high_exponent_real<D, E, N1> a,
                                          high_exponent_real<D, E, N2> b) {
  return a + -b;
}

template <typename D, typename E, bool N1, bool N2>
high_exponent_real<D, E, false> operator+(high_exponent_real<D, E, N1> a,
                                          high_exponent_real<D, E, N2> b) {
  if (a.mantissa() == 0)
    return b;
  if (b.mantissa() == 0)
    return a;

  // Need to figure out a sensible exponent
  if (a.exponent() > b.exponent()) {
    // Convert b to a's exponent
    return {a.mantissa() +
                b.mantissa() * (D)std::exp2(b.exponent() - a.exponent()),
            a.exponent()};
  } else if (a.exponent() < b.exponent()) {
    // Convert a to b's exponent
    return {a.mantissa() * (D)std::exp2(a.exponent() - b.exponent()) +
                b.mantissa(),
            b.exponent()};
  } else {
    return {a.mantissa() + b.mantissa(), a.exponent()};
  }
}

template <typename D, typename E, bool N1, bool N2>
high_exponent_real<D, E, false> operator*(high_exponent_real<D, E, N1> a,
                                          high_exponent_real<D, E, N2> b) {
  return {a.mantissa() * b.mantissa(), a.exponent() + b.exponent()};
}

template <typename D, typename E, bool N1, bool N2>
high_exponent_real<D, E, false> operator/(high_exponent_real<D, E, N1> a,
                                          high_exponent_real<D, E, N2> b) {
  return {a.mantissa() / b.mantissa(), a.exponent() - b.exponent()};
}

template <typename D, typename E>
std::ostream &operator<<(std::ostream &os, high_exponent_real<D, E> a) {
  os << a.mantissa();
  if (a.exponent())
    os << "*2^" << a.exponent();
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
  if (a.mantissa() < 0 && b.mantissa() > 0)
    return -1;
  if (a.mantissa() > 0 && b.mantissa() < 0)
    return 1;

  int a_bigger = a.mantissa() < 0 ? -1 : 1;

  // !! This does not work for zero !! FIXME
  // They are the same sign
  // Compare exponents
  if (a.exponent() > b.exponent())
    return a_bigger;
  if (a.exponent() < b.exponent())
    return -a_bigger;

  // The exponent is the same
  if (a.mantissa() < b.mantissa())
    return -1;
  if (a.mantissa() > b.mantissa())
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
  return std::log(a.mantissa()) + a.exponent() * std::log(2);
}

template <typename D, typename E, bool N>
high_exponent_real<D, E, false> operator<<(const high_exponent_real<D, E, N> &a,
                                           int shift) {
  return {a.mantissa(), a.exponent() + shift};
}

template <typename D, typename E, bool N>
high_exponent_real<D, E, false> operator>>(const high_exponent_real<D, E, N> &a,
                                           int shift) {
  return a << -shift;
}

// !! Move into complex_number.hpp
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

using high_exponent_double = high_exponent_real<double, int, false>;

} // namespace fractals
