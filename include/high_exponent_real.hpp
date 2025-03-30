#pragma once

#include <cmath>
#include <iostream>

namespace fractals {

/*
  A high_exponent_real is a normal floating point number (like a double) with
  extra precision for the exponent.
*/
template <std::floating_point Double = double, std::integral Exponent = int,
          bool Normalized = false>
class high_exponent_real {
public:
  template <bool N>
  constexpr high_exponent_real(
      const high_exponent_real<Double, Exponent, N> &src) {
    *this = src;
  }

  template <bool N>
  high_exponent_real &
  operator=(const high_exponent_real<Double, Exponent, N> &src) {
    if constexpr (!N) {
      assign(src.mantissa(), src.exponent());
    } else {
      d = src.mantissa();
      e = src.exponent();
    }
    return *this;
  }

  high_exponent_real<Double, Exponent, true> normalize() const { return *this; }

  constexpr void assign(Double d0, Exponent e0) {
    if constexpr (Normalized) {
      d = std::frexp(d0, &e);
      e += e0;
    } else {
      d = d0;
      e = e0;
    }
  }

  constexpr high_exponent_real(Double d0 = 0, Exponent e0 = 0) {
    assign(d0, e0);
  }

  Double to_double() const { return std::ldexp(d, e); }
  Double mantissa() const { return d; }
  Exponent exponent() const { return e; }

  template <bool N>
  high_exponent_real &operator+=(high_exponent_real<Double, Exponent, N> a) {
    return *this = (*this + a);
  }

  template <bool N>
  high_exponent_real &operator-=(high_exponent_real<Double, Exponent, N> a) {
    return *this = (*this - a);
  }

  template <bool N>
  high_exponent_real &operator*=(high_exponent_real<Double, Exponent, N> a) {
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
high_exponent_real<D, E, N> operator-(high_exponent_real<D, E, N> a) {
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

template <typename D, typename E, bool N>
high_exponent_real<D, E, false> operator*(high_exponent_real<D, E, N> a,
                                          int b) {
  return {a.mantissa() * b, a.exponent()};
}

template <typename D, typename E, bool N>
high_exponent_real<D, E, false> operator*(int a,
                                          high_exponent_real<D, E, N> b) {
  return {a * b.mantissa(), b.exponent()};
}

template <typename D, typename E, bool N1, bool N2>
high_exponent_real<D, E, false> operator/(high_exponent_real<D, E, N1> a,
                                          high_exponent_real<D, E, N2> b) {
  return {a.mantissa() / b.mantissa(), a.exponent() - b.exponent()};
}

template <typename D, typename E, bool N>
std::ostream &operator<<(std::ostream &os, high_exponent_real<D, E, N> a) {
  os << a.mantissa();
  if (a.exponent())
    os << "*2^" << a.exponent();
  return os;
}

namespace detail {
template <typename D, typename E>
int cmp(high_exponent_real<D, E, true> a, high_exponent_real<D, E, true> b) {
  // If the sign is different
  if (a.mantissa() <= 0 && b.mantissa() > 0)
    return -1;
  if (a.mantissa() >= 0 && b.mantissa() < 0)
    return 1;

  int a_bigger = a.mantissa() < 0 ? -1 : 1;

  if (b.mantissa() == 0) {
    if (a.mantissa() == 0)
      return 0;
    return a_bigger;
  }

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
} // namespace detail

template <typename D, typename E, bool N1, bool N2>
bool operator==(high_exponent_real<D, E, N1> a,
                high_exponent_real<D, E, N2> b) {
  return detail::cmp(normalize(a), normalize(b)) == 0;
}

template <typename D, typename E, bool N1, bool N2>
bool operator>=(high_exponent_real<D, E, N1> a,
                high_exponent_real<D, E, N2> b) {
  return detail::cmp(normalize(a), normalize(b)) >= 0;
}

template <typename D, typename E, bool N1, bool N2>
bool operator>(high_exponent_real<D, E, N1> a, high_exponent_real<D, E, N2> b) {
  return detail::cmp(normalize(a), normalize(b)) > 0;
}

template <typename D, typename E, bool N1, bool N2>
bool operator<(high_exponent_real<D, E, N1> a, high_exponent_real<D, E, N2> b) {
  return detail::cmp(normalize(a), normalize(b)) < 0;
}

template <typename D, typename E, bool N1, bool N2>
bool operator<=(high_exponent_real<D, E, N1> a,
                high_exponent_real<D, E, N2> b) {
  return detail::cmp(normalize(a), normalize(b)) <= 0;
}

template <typename D, typename E, bool N>
bool isfinite(high_exponent_real<D, E, N> a) {
  return std::isfinite(a.mantissa());
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

using high_exponent_double = high_exponent_real<double, int, false>;

} // namespace fractals
