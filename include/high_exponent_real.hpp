#pragma once
#include "convert.hpp"
#include "high_precision_real.hpp"
#include "real_number.hpp"

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
};

template <typename D, typename E>
high_exponent_real<D, E> normalize(const high_exponent_real<D, E> e) {
  return {e.d, e.e};
}

template <typename D, typename E>
high_exponent_real<D, E> operator-(high_exponent_real<D, E> a) {
  return {-a.d, a.e, skip_normalization()};
}

template <typename D, typename E>
high_exponent_real<D, E> operator-(high_exponent_real<D, E> a,
                                   high_exponent_real<D, E> b) {
  return a + -b;
}

template <typename D, typename E>
high_exponent_real<D, E> operator+(high_exponent_real<D, E> a,
                                   high_exponent_real<D, E> b) {
  if (a.d == 0)
    return b;
  if (b.d == 0)
    return a;

  // Need to figure out a sensible exponent
  if (a.e > b.e) {
    // Convert b to a's exponent
    return {a.d + b.d * (D)std::exp2(b.e - a.e), a.e, skip_normalization()};
  } else if (a.e < b.e) {
    // Convert a to b's exponent
    return {a.d * (D)std::exp2(a.e - b.e) + b.d, b.e, skip_normalization()};
  } else {
    return {a.d + b.d, a.e, skip_normalization()};
  }
}

template <typename D, typename E>
high_exponent_real<D, E> operator*(high_exponent_real<D, E> a,
                                   high_exponent_real<D, E> b) {
  return {a.d * b.d, a.e + b.e, skip_normalization()};
}

template <typename D, typename E>
high_exponent_real<D, E> operator/(high_exponent_real<D, E> a,
                                   high_exponent_real<D, E> b) {
  return high_exponent_real<D, E>{a.d / b.d, a.e - b.e, skip_normalization()};
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
int cmp(high_exponent_real<D, E> a, high_exponent_real<D, E> b) {
  a = normalize(a);
  b = normalize(b);

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
  a = normalize(a);
  b = normalize(b);
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

template <std::floating_point D1, std::floating_point D2, typename E>
struct convert_to<D1, high_exponent_real<D2, E>> {
  static D1 get(const high_exponent_real<D2, E> &x) { return x.to_double(); }
};

template <std::floating_point D1, std::floating_point D2, typename E>
struct convert_to<high_exponent_real<D1, E>, D2> {
  static high_exponent_real<D1, E> get(D2 x) { return {(D1)x}; }
};

template <typename D, typename E>
struct convert_to<high_exponent_real<D, E>, int> {
  static high_exponent_real<D, E> get(int x) { return {x}; }
};

template <int N, typename D, typename E>
struct convert_to<high_precision_real<N>, high_exponent_real<D, E>> {
  static high_precision_real<N> get(const high_exponent_real<D, E> &x) {
    return high_precision_real<N>(x.d) << x.e;
  }
};

template <int N, std::floating_point D, typename E>
struct convert_to<high_exponent_real<D, E>, high_precision_real<N>> {
  static high_exponent_real<D, E> get(const high_precision_real<N> &x) {
    if (x.fraction[0] == 0) {
      int e = count_fractional_zeros(x);
      return {(D)(x << e).to_double(), -e};
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

template <int Digits, int MinExp, int MaxExp>
  requires(Digits <= std::numeric_limits<long double>::digits &&
           (MinExp<std::numeric_limits<long double>::min_exponent || MaxExp>
                std::numeric_limits<long double>::max_exponent))
struct make_real<Digits, MinExp, MaxExp> {
  using type = high_exponent_real<long double, int>;
};

// Specialise this because std::complex does incompatible things
template <typename D, typename E>
std::complex<high_exponent_real<D, E>>
operator*(std::complex<high_exponent_real<D, E>> a,
          std::complex<high_exponent_real<D, E>> b) {
  return {a.real() * b.real() - a.imag() * b.imag(),
          a.real() * b.imag() + a.imag() * b.real()};
}

} // namespace fractals
