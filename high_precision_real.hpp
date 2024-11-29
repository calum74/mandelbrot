#pragma once

#include <cmath>
#include <cstdint>
#include <iostream>

#define HP_FLOAT_VALIDATION 0

#if HP_FLOAT_VALIDATION
#include <cassert>
#endif
#include "convert.hpp"

namespace fractals {

/*
    A minimal high-precision float implementation for greater precision.
*/
template <int N> struct high_precision_real {

  static_assert(N > 0);

  // The sign of the number
  bool negative;

  // Fixed-point binary, with fraction[0] representing the integer component,
  // and fraction[1..] representing the fractional parts to the given precision
  std::uint64_t fraction[N];

  high_precision_real() : negative(false), fraction{} {}

  high_precision_real(int n) : negative(n < 0), fraction{} {
    fraction[0] = negative ? -n : n;
  }

  template <int M>
  high_precision_real(const high_precision_real<M> &other)
      : negative(other.negative) {
    for (int i = 0; i < N; ++i) {
      fraction[i] = i >= M ? 0 : other.fraction[i];
    }
  }

  high_precision_real &operator+=(const high_precision_real &b) {
    // This could be implemented more efficiently
    return *this = *this + b;
  }

  high_precision_real &operator*=(const high_precision_real &b) {
    // This could be implemented more efficiently
    return *this = *this * b;
  }

  high_precision_real(double d) : negative(d < 0) {
    if (negative)
      d = -d;
    for (int i = 0; i < N; ++i) {
      double ip;
      d = std::modf(d, &ip);
      fraction[i] = ip;
      d *= std::pow(2.0, 64);
    }
  }

  double to_double() const {
    double value = fraction[0];
    for (int i = 1; i < N; ++i)
      value += fraction[i] * pow(0.5, 64 * i);
    return negative ? -value : value;
  }

  int to_int() const { return negative ? -fraction[0] : fraction[0]; }
};

template <int N>
void raw_lshift(const high_precision_real<N> &a,
                high_precision_real<N> &result) {
  int carry = 0;
  for (int i = N - 1; i >= 0; --i) {
    int new_carry = (a.fraction[i] & (1ull << 63)) != 0;
    result.fraction[i] = (a.fraction[i] << 1) | carry;
    carry = new_carry;
  }
}

template <int N>
void raw_rshift(const high_precision_real<N> &a,
                high_precision_real<N> &result) {
  std::uint64_t carry = 0;
  for (int i = 0; i < N; ++i) {
    std::uint64_t new_carry = (a.fraction[i] & 1) ? (1ull << 63) : 0;
    result.fraction[i] = (a.fraction[i] >> 1) | carry;
    carry = new_carry;
  }
}

template <int N> double to_double(const high_precision_real<N> &a) {
  return a.to_double();
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

template <int M> struct convert_to<high_precision_real<M>, double> {
  static high_precision_real<M> get(double x) { return {x}; }
};

template <int N>
high_precision_real<N> operator/(const high_precision_real<N> &a, double d) {
  return a * high_precision_real<N>{1.0 / d};
}

template <int N>
high_precision_real<N> operator/(const high_precision_real<N> &a, int n) {
  if (n == 2) {
    high_precision_real<N> result;
    raw_rshift(a, result);
    result.negative = a.negative;
    return result;
  } else
    return a * high_precision_real<N>{1.0 / n};
}

template <int N>
int cmp(const high_precision_real<N> &a, const high_precision_real<N> &b) {
  if (a.negative && !b.negative)
    return -1;
  if (!a.negative && b.negative)
    return 1;
  return a.negative ? -raw_cmp(a, b) : raw_cmp(a, b);
}

template <int N>
bool operator<(const high_precision_real<N> &a,
               const high_precision_real<N> &b) {
  return cmp(a, b) < 0;
}

template <int N> bool operator>=(const high_precision_real<N> &a, int b) {
  // TODO: Doesn't handle negative cases properly
  if (a.negative)
    return false;
  return a.fraction[0] >= b;
}

template <int N> bool operator<(const high_precision_real<N> &a, int b) {
  // TODO: Doesn't handle negative cases properly
  if (a.negative)
    return false;
  return a.fraction[0] < b;
}

template <int N> bool operator<=(const high_precision_real<N> &a, int b) {
  // TODO: Doesn't handle negative cases properly
  if (a.negative)
    return false;
  return a.fraction[0] <= b;
}

template <int N> high_precision_real<N> operator-(high_precision_real<N> b) {
  b.negative = !b.negative;
  return b;
}

template <int N>
high_precision_real<N> inverse(const high_precision_real<N> &b) {

  // TODO!!!
  return b;
}

template <int N>
bool operator==(const high_precision_real<N> &a,
                const high_precision_real<N> &b) {
  return cmp(a, b) == 0;
}

template <int N>
high_precision_real<N> operator*(const high_precision_real<N> &a,
                                 const high_precision_real<N> &b) {
  high_precision_real<N> result;
  raw_mul(a, b, result);
  result.negative = a.negative != b.negative;
#if HP_FLOAT_VALIDATION
  auto a_debug = a.to_double();
  auto b_debug = b.to_double();
  auto c_debug = result.to_double();
  assert(c_debug == 0 || abs(a_debug * b_debug / c_debug - 1) < 0.001);
#endif
  return result;
}

template <int N>
high_precision_real<N> operator*(const high_precision_real<N> &a, double b) {
  return a * high_precision_real<N>{b};
}

template <int N>
high_precision_real<N> operator*(int a, const high_precision_real<N> &b) {
  if (a == 2) {
    high_precision_real<N> result;
    raw_lshift(b, result);
    result.negative = b.negative;
    return result;
  } else {
    return high_precision_real<N>{a} * b;
  }
}

template <int N>
void raw_add(const high_precision_real<N> &a, const high_precision_real<N> &b,
             high_precision_real<N> &result) {
  result.fraction[N - 1] = a.fraction[N - 1] + b.fraction[N - 1];

  for (int i = N - 2; i >= 0; i--) {
    result.fraction[i] = a.fraction[i] + b.fraction[i] +
                         (result.fraction[i + 1] < b.fraction[i + 1]);
  }
}

template <int N>
void raw_mul(const high_precision_real<N> &a, const high_precision_real<N> &b,
             high_precision_real<N> &result) {
  for (int i = N - 1; i >= 0; i--)
    for (int j = N - 1 - i; j >= 0; j--) {
      auto m = (__uint128_t)a.fraction[i] * (__uint128_t)b.fraction[j];
      std::uint64_t m1 = m;
      std::uint64_t m2 = m >> 64;
      auto ij = i + j;
      result.fraction[ij] += m1;
      int carry = result.fraction[ij] < m1;
      m2 += carry;
      if (ij > 0) {
        ij--;
        result.fraction[ij] += m2;
        carry = result.fraction[ij] < m2;
      }
      while (carry && ij-- > 0) {
        result.fraction[ij] += carry;
        carry = result.fraction[ij] == 0;
      }
    }
}

template <int N>
int raw_cmp(const high_precision_real<N> &a, const high_precision_real<N> &b) {
  for (int i = 0; i < N; ++i) {
    if (a.fraction[i] < b.fraction[i])
      return -1;
    if (a.fraction[i] > b.fraction[i])
      return 1;
  }
  return 0;
}

template <int N>
void raw_sub(const high_precision_real<N> &a, const high_precision_real<N> &b,
             high_precision_real<N> &result) {
  result.fraction[N - 1] = a.fraction[N - 1] - b.fraction[N - 1];

  for (int i = N - 2; i >= 0; i--) {
    result.fraction[i] = a.fraction[i] - b.fraction[i] -
                         (result.fraction[i + 1] > a.fraction[i + 1]);
  }
}

template <int N>
high_precision_real<N> operator+(const high_precision_real<N> &a,
                                 const high_precision_real<N> &b) {
  high_precision_real<N> result;

  if (a.negative == b.negative) {
    result.negative = a.negative;
    raw_add(a, b, result);
  } else {
    int c = raw_cmp(a, b);
    if (c < 0) {
      // a<b:
      raw_sub(b, a, result);
      result.negative = b.negative;
    } else {
      // a>b
      raw_sub(a, b, result);
      result.negative = a.negative;
    }
  }

#if HP_FLOAT_VALIDATION
  auto a_debug = a.to_double();
  auto b_debug = b.to_double();
  auto r_debug = result.to_double();
  assert(abs(a_debug + b_debug - r_debug) < 0.001);
#endif
  return result;
}

template <int N>
high_precision_real<N> operator-(const high_precision_real<N> &a,
                                 const high_precision_real<N> &b) {
  high_precision_real<N> result;

  if (a.negative != b.negative) {
    result.negative = a.negative;
    raw_add(a, b, result);
  } else {
    int c = raw_cmp(a, b);
    if (c < 0) {
      // a<b:
      raw_sub(b, a, result);
      result.negative = !a.negative;
    } else {
      // a>b
      raw_sub(a, b, result);
      result.negative = a.negative;
    }
  }

#if HP_FLOAT_VALIDATION
  auto a_debug = a.to_double();
  auto b_debug = b.to_double();
  auto r_debug = result.to_double();
  assert(abs(a_debug - b_debug - r_debug) < 0.001);
#endif

  return result;
}

template <int N>
std::ostream &operator<<(std::ostream &os, high_precision_real<N> n) {

  if (n.negative)
    os << '-';

  if (os.flags() & std::ios_base::hex) {
    os << "0x";
    // TODO: Hex mode
  }

  os << n.fraction[0] << '.';

  if (os.flags() & std::ios_base::hex) {
    // TODO: Hex output
  } else {
    int digits = (N - 1) * 64 * 0.301; // 0.301 = std::log(2) / std::log(10);

    // A naive O(N^2) algorithm.
    high_precision_real<N> ten{10};

    if (digits > os.precision())
      digits = os.precision();

    for (int i = 0; i < digits; i++) {
      n.fraction[0] = 0;
      n = n * ten;
      os << n.fraction[0];
    }
  }

  return os;
}

template <int N>
std::istream &operator>>(std::istream &is, high_precision_real<N> &n) {

  char ch;
  bool negative;
  ch = is.peek();
  n.negative = false;
  if (ch == '-') {
    is.get(ch);
    negative = true;
  } else {
    negative = false;
  }

  // Read up to the decimal point

  for (int i = 0; i < N; i++)
    n.fraction[i] = 0;

  do {
    ch = is.peek();
    if (std::isdigit(ch)) {
      is.get(ch);
      n.fraction[0] = n.fraction[0] * 10 + (ch - '0');
    }
  } while (std::isdigit(ch));

  if (ch != '.') {
    n.negative = negative;
    return is;
  }
  is.get(ch); // Skip '.'

  // Construct 0.1 in binary
  high_precision_real<N> tenth, mult;
  tenth.fraction[1] = 0x1999999999999999ull;
  for (int i = 2; i < N; ++i) {
    tenth.fraction[i] = 0x9999999999999999ull;
  }
  mult = tenth;

  do {
    ch = is.peek();
    if (std::isdigit(ch)) {
      is.get(ch);
      high_precision_real<N> i;
      // Not very optimal algorithm
      i.fraction[0] = (ch - '0');
      n = n + (i * mult);
      mult = mult * tenth;
    }
  } while (std::isdigit(ch));
  n.negative = negative; // Set last because we're adding in the loop
  return is;
}

template <int N> int count_zeros(const high_precision_real<N> &n) {
  int c = 0;
  for (int i = 1; i < N; ++i, c += 64) {
    if (n.fraction[i]) {
      for (std::uint64_t b = (std::uint64_t)1 << 63; b; c++, b >>= 1) {
        if (b & n.fraction[i]) {
          return c;
        }
      }
    }
  }
  return 0;
}

} // namespace fractals