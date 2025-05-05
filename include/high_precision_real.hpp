/*
  Naive implementation of high-precision arithmetic in C++.

  These algorithms are not particularly efficient, for example multiplication is
  quadratic.

  It's also not very well tested, so don't use this for anything other than
  displaying pretty pictures.
*/

#pragma once

#include "number_cast.hpp"

#include <cmath>
#include <cstdint>
#include <iostream>

#if _WIN32
#include <intrin.h>
#if defined _M_ARM64
#pragma intrinsic(__umulh)
#else
#pragma intrinsic(_umul128)
#endif
#endif

namespace numbers
{

  /*
      A simple high-precision number implementation.
      The number of integer bits is 64 (part_type), and the number of fractional
      bits is FractionalBits.
  */
  template <int FractionalBits>
    requires(FractionalBits >= 0)
  class high_precision_real
  {
  public:
    using part_type = std::uint64_t;
    static constexpr int part_bits = sizeof(part_type) * 8; // 8 bits per byte
    static constexpr int parts = 1 + (FractionalBits + part_bits - 1) / part_bits;
    using size_type = int;

    // Initialize to 0
    high_precision_real() : sign(false), part{} {}

    // Initialize from an integer
    high_precision_real(int n) : sign(n < 0), part{}
    {
      part[0] = negative() ? -n : n;
    }

    // Initialize from another high_precision_real, truncates or pads the extra
    // bits
    template <int M>
    high_precision_real(const high_precision_real<M> &other)
        : sign(other.negative())
    {
      for (int i = 0; i < size(); ++i)
      {
        part[i] = i >= other.size() ? 0 : other[i];
      }
    }

    high_precision_real &operator*=(const high_precision_real &b)
    {
      // This could be implemented more efficiently
      return *this = *this * b;
    }

    high_precision_real(double d) : sign(d < 0)
    {
      if (negative())
        d = -d;
      for (int i = 0; i < parts; ++i)
      {
        double ip;
        d = std::modf(d, &ip);
        part[i] = ip;
        d *= std::pow(2.0, part_bits);
      }
    }

    high_precision_real &operator+=(const high_precision_real &b)
    {
      // This could be implemented more efficiently
      return *this = *this + b;
    }

    high_precision_real &operator-=(const high_precision_real &b)
    {
      // This could be implemented more efficiently
      return *this = *this - b;
    }

    double to_double() const
    {
      double value = part[0];
      for (int i = 1; i < size(); ++i)
        value += part[i] * std::pow(0.5, part_bits * i);
      return negative() ? -value : value;
    }

    // Access to the parts
    // [0] is the integer component, [1:] are the fractional parts
    part_type &operator[](int i) { return part[i]; }
    const part_type &operator[](int i) const { return part[i]; }

    // The total number of parts
    constexpr size_type size() const { return parts; }

    // Get/set if the number is negative
    bool negative() const { return sign; }
    bool &negative() { return sign; }

  private:
    // The sign of the number: 0 = positive, 1 = negative
    bool sign;

    // Fixed-point binary, with part[0] representing the integer component,
    // and part[1:] representing the fractional parts
    part_type part[parts];
  };

  namespace detail
  {
    template <int N>
    void raw_lshift(const high_precision_real<N> &a,
                    high_precision_real<N> &result)
    {
      int carry = 0;
      constexpr int part_bits = high_precision_real<N>::part_bits;
      for (int i = a.size() - 1; i >= 0; --i)
      {
        int new_carry = (a[i] & (1ull << (part_bits - 1))) != 0;
        result[i] = (a[i] << 1) | carry;
        carry = new_carry;
      }
    }

    template <int N>
    void raw_rshift(const high_precision_real<N> &a,
                    high_precision_real<N> &result)
    {
      using part_type = typename high_precision_real<N>::part_type;
      constexpr int part_bits = high_precision_real<N>::part_bits;
      part_type carry = 0;
      for (int i = 0; i < a.size(); ++i)
      {
        part_type new_carry = (a[i] & 1) ? (1ull << (part_bits - 1)) : 0;
        result[i] = (a[i] >> 1) | carry;
        carry = new_carry;
      }
    }

    template <int N>
    int raw_cmp(const high_precision_real<N> &a, const high_precision_real<N> &b);

    template <int N>
    int cmp(const high_precision_real<N> &a, const high_precision_real<N> &b)
    {
      if (is_zero(a) && is_zero(b))
        return 0;
      if (a.negative() && !b.negative())
        return -1;
      if (!a.negative() && b.negative())
        return 1;
      return a.negative() ? -raw_cmp(a, b) : raw_cmp(a, b);
    }

    template <int N>
    void raw_add(const high_precision_real<N> &a, const high_precision_real<N> &b,
                 high_precision_real<N> &result)
    {
      result[result.size() - 1] = a[a.size() - 1] + b[b.size() - 1];

      for (int i = a.size() - 2; i >= 0; i--)
      {
        result[i] = a[i] + b[i] + (result[i + 1] < a[i + 1]);
      }
    }

    template <int N>
    void raw_mul(const high_precision_real<N> &a, const high_precision_real<N> &b,
                 high_precision_real<N> &result)
    {
      for (int i = a.size() - 1; i >= 0; i--)
      {
        // TODO: Fix bounds of inner loop
        for (int j = a.size() - 1; j >= 0; j--)
        {
#if _WIN32
#if defined _M_ARM64
          std::uint64_t m1 = a[i] * b[j];
          std::uint64_t m2 = __umulh(a[i], b[j]);
#else
          std::uint64_t m2;
          std::uint64_t m1 = _umul128(a[i], b[j], &m2);
#endif
#else
          __uint128_t m = (__uint128_t)a[i] * (__uint128_t)b[j];
          std::uint64_t m1 = m;
          std::uint64_t m2 = m >> 64;
#endif
          auto ij = i + j;
          int carry;
          if (ij < a.size())
          {
            result[ij] += m1;
            carry = result[ij] < m1;
            m2 += carry;
            carry = m2 < carry;
          }
          else
            carry = 0;
          if (ij > 0 && ij <= a.size())
          {
            ij--;
            result[ij] += m2;
            carry = result[ij] < m2;
          }
          // Move into loop TODO
          while (carry && ij-- > 0)
          {
            result[ij] += carry;
            carry = result[ij] == 0;
          }
        }
      }
    }

    template <int N>
    int raw_cmp(const high_precision_real<N> &a, const high_precision_real<N> &b)
    {
      for (int i = 0; i < a.size(); ++i)
      {
        if (a[i] < b[i])
          return -1;
        if (a[i] > b[i])
          return 1;
      }
      return 0;
    }

    template <int N>
    void raw_sub(const high_precision_real<N> &a, const high_precision_real<N> &b,
                 high_precision_real<N> &result)
    {
      std::uint64_t borrow = 0;

      for (int i = a.size() - 1; i >= 0; i--)
      {
        result[i] = a[i] - b[i] - borrow;
        borrow = a[i] + borrow < b[i];
      }
    }

    template <int N>
    void raw_shiftleft(const high_precision_real<N> &a, high_precision_real<N> &b,
                       int n)
    {
      using part_type = typename high_precision_real<N>::part_type;
      constexpr int part_bits = high_precision_real<N>::part_bits;
      part_type extra = 0;
      int m = n / part_bits;
      n = n % part_bits;
      for (int i = a.size() - 1; i >= 0; i--)
      {
        if (i + m < a.size())
        {
          b[i] = (a[i + m] << n) | extra;
          extra = n == 0 ? 0 : a[i + m] >> (part_bits - n);
        }
        else
          b[i] = 0;
      }
    }

    template <int N>
    void raw_shiftright(const high_precision_real<N> &a, high_precision_real<N> &b,
                        int n)
    {
      using part_type = typename high_precision_real<N>::part_type;
      part_type extra = 0;
      constexpr int part_bits = high_precision_real<N>::part_bits;
      int m = n / part_bits;
      n = n % part_bits;
      for (int i = 0; i < a.size(); i++)
      {
        if (i - m >= 0)
        {
          b[i] = a[i - m] >> n | extra;
          extra = n == 0 ? 0 : a[i - m] << (part_bits - n);
        }
        else
          b[i] = 0;
      }
    }
  } // namespace detail

  template <int N>
  high_precision_real<N> operator/(const high_precision_real<N> &a, double d)
  {
    return a * high_precision_real<N>{1.0 / d};
  }

  template <int N>
  high_precision_real<N> operator/(const high_precision_real<N> &a, int n)
  {
    if (n == 2)
    {
      high_precision_real<N> result;
      detail::raw_rshift(a, result);
      result.negative() = a.negative();
      return result;
    }
    else
      return a * high_precision_real<N>{1.0 / n};
  }

  template <int N>
  bool operator<(const high_precision_real<N> &a,
                 const high_precision_real<N> &b)
  {
    return detail::cmp(a, b) < 0;
  }

  template <int N>
  bool operator>=(const high_precision_real<N> &a,
                  const high_precision_real<N> &b)
  {
    return detail::cmp(a, b) >= 0;
  }

  template <int N>
  bool operator>=(const high_precision_real<N> &a, int b)
  {
    // TODO: Doesn't handle negative cases properly
    if (a.negative())
      return false;
    return a[0] >= b;
  }

  template <int N>
  bool operator<(const high_precision_real<N> &a, int b)
  {
    // TODO: Doesn't handle negative cases properly
    if (a.negative)
      return false;
    return a.fraction[0] < b;
  }

  template <int N>
  bool operator<=(const high_precision_real<N> &a, int b)
  {
    // TODO: Doesn't handle negative cases properly
    if (a.negative())
      return true;
    if (a[0] < b)
      return true;
    if (a[0] > b)
      return false;
    for (int i = 1; i < a.size(); i++)
      if (a[i])
        return false;
    return true;
  }

  template <int N>
  high_precision_real<N> operator-(high_precision_real<N> b)
  {
    if (!is_zero(b))
      b.negative() = !b.negative();
    return b;
  }

  template <int N>
  int count_leading_zeros(const high_precision_real<N> &n)
  {
    int c = 0;
    for (int i = 0; i < n.size(); ++i, c += 64)
    {
      if (n[i])
      {
        for (std::uint64_t b = (std::uint64_t)1 << 63; b; c++, b >>= 1)
        {
          if (b & n[i])
          {
            return c;
          }
        }
      }
    }
    return c;
  }

  template <int N>
  high_precision_real<N> inverse(const high_precision_real<N> &d)
  {

    // See
    // https://en.wikipedia.org/wiki/Division_algorithm#Newton%E2%80%93Raphson_division

    // Step 1: Scale d such that it's beween 0.5 and less than 1
    // 1a, find out the number of leading zeros of b
    int m = count_leading_zeros(d) - 64;
    auto D2 = d << m; // Shift d left or right
    D2.negative() = false;

    // Step 2: Compute X = 48/17 − 32/17 × D'
    high_precision_real<N> X = high_precision_real<N>{48.0 / 17.0} -
                               D2 * high_precision_real<N>{32.0 / 17.0};
    int iterations = std::ceil(std::log2(64 * N / std::log2(17)));
    high_precision_real<N> one = 1;
    for (int i = 0; i < iterations; ++i)
    {
      X = X + X * (one - D2 * X);
    }
    auto result = X << m;
    result.negative() = d.negative();
    return result;
  }

  template <int N>
  bool operator==(const high_precision_real<N> &a,
                  const high_precision_real<N> &b)
  {
    return detail::cmp(a, b) == 0;
  }

  template <int N>
  bool is_zero(const high_precision_real<N> &a)
  {
    for (int i = 0; i < a.size(); ++i)
      if (a[i])
        return false;
    return true;
  }

  template <int N>
  high_precision_real<N> operator*(const high_precision_real<N> &a, double b)
  {
    return a * high_precision_real<N>{b};
  }

  template <int N>
  high_precision_real<N> operator*(int a, const high_precision_real<N> &b)
  {
    if (a == 2)
    {
      high_precision_real<N> result;
      detail::raw_lshift(b, result);
      result.negative() = b.negative();
      return result;
    }
    else
    {
      return high_precision_real<N>{a} * b;
    }
  }

  template <int N>
  high_precision_real<N> operator*(const high_precision_real<N> &a,
                                   const high_precision_real<N> &b)
  {
    high_precision_real<N> result;
    detail::raw_mul(a, b, result);
    result.negative() = a.negative() != b.negative();
    if (is_zero(result))
      result.negative() = false;
    return result;
  }

  template <int N>
  high_precision_real<N> operator+(const high_precision_real<N> &a,
                                   const high_precision_real<N> &b)
  {
    high_precision_real<N> result;

    if (a.negative() == b.negative())
    {
      result.negative() = a.negative();
      detail::raw_add(a, b, result);
    }
    else
    {
      int c = detail::raw_cmp(a, b);
      if (c < 0)
      {
        // a<b:
        detail::raw_sub(b, a, result);
        result.negative() = b.negative();
      }
      else if (c > 0)
      {
        // a>b
        detail::raw_sub(a, b, result);
        result.negative() = a.negative();
      }
      // else: result = 0
    }

    return result;
  }

  template <int N>
  high_precision_real<N> operator-(const high_precision_real<N> &a,
                                   const high_precision_real<N> &b)
  {
    high_precision_real<N> result;

    if (a.negative() != b.negative())
    {
      result.negative() = a.negative();
      detail::raw_add(a, b, result);
    }
    else
    {
      int c = detail::raw_cmp(a, b);
      if (c < 0)
      {
        // a<b:
        detail::raw_sub(b, a, result);
        result.negative() = !a.negative();
      }
      else if (c > 0)
      {
        // a>b
        detail::raw_sub(a, b, result);
        result.negative() = a.negative();
      }
      // else: 0
    }

    return result;
  }

  template <int N>
  std::ostream &operator<<(std::ostream &os, high_precision_real<N> n)
  {

    if (n.negative())
      os << '-';

    if (os.flags() & std::ios_base::hex)
    {
      os << "0x";
      // TODO: Hex mode
    }

    os << n[0] << '.';

    if (os.flags() & std::ios_base::hex)
    {
      // TODO: Hex output
      os << "<todo>";
    }
    else
    {
      int digits = (n.size() - 1) * high_precision_real<N>::part_bits *
                   0.301; // 0.301 = std::log(2) / std::log(10);

      // A naive O(N^2) algorithm.
      high_precision_real<N> ten{10};

      if (digits > os.precision())
        digits = os.precision();

      for (int i = 0; i < digits; i++)
      {
        n[0] = 0;
        n = n * ten;
        os << n[0];
      }
    }

    return os;
  }

  template <int N>
  void make_tenth(high_precision_real<N> &tenth)
  {
    tenth.negative() = false;
    tenth[0] = 0;
    tenth[1] = 0x1999999999999999ull;
    for (int i = 2; i < tenth.size() - 1; ++i)
    {
      tenth[i] = 0x9999999999999999ull;
    }
    if (tenth.size() > 2)
      tenth[tenth.size() - 1] = 0x999999999999999aull;
  }

  template <int N>
  std::istream &operator>>(std::istream &is, high_precision_real<N> &n)
  {
    char ch;
    bool negative;
    ch = is.peek();
    n.negative() = false;
    if (ch == '-')
    {
      is.get(ch);
      negative = true;
    }
    else
    {
      negative = false;
    }

    // Read up to the decimal point

    for (int i = 0; i < n.size(); i++)
      n[i] = 0;

    do
    {
      ch = is.peek();
      if (std::isdigit(ch))
      {
        is.get(ch);
        n[0] = n[0] * 10 + (ch - '0');
      }
    } while (std::isdigit(ch));

    // Construct 0.1 in binary
    high_precision_real<N> tenth, mult;
    make_tenth(tenth);
    mult = tenth;

    if (ch != 'e')
    {
      if (ch != '.')
      {
        n.negative() = negative;
        return is;
      }
      is.get(ch); // Skip '.'

      do
      {
        ch = is.peek();
        if (std::isdigit(ch))
        {
          is.get(ch);
          high_precision_real<N> i;
          // Not very optimal algorithm
          i[0] = (ch - '0');
          n = n + (i * mult);
          mult = mult * tenth;
        }
      } while (std::isdigit(ch));
    }

    if (ch == 'e')
    {
      int power = 0;
      is >> ch >> power;

      // This is super-inefficient
      while (power < 0)
      {
        n = n * tenth;
        power++;
      }
      while (power > 0)
      {
        n = n * 10;
        power--;
      }
    }

    n.negative() = negative; // Set last because we're adding in the loop
    return is;
  }

  template <int N>
  int count_fractional_zeros(const high_precision_real<N> &n)
  {
    int c = 0;
    using part_type = typename high_precision_real<N>::part_type;
    const int part_bits = high_precision_real<N>::part_bits;
    for (int i = 1; i < n.size(); ++i, c += part_bits)
    {
      if (n[i])
      {
        for (part_type b = (part_type)1 << (part_bits - 1); b; c++, b >>= 1)
        {
          if (b & n[i])
          {
            return c;
          }
        }
      }
    }
    return 0;
  }

  template <int N>
  high_precision_real<N> operator<<(const high_precision_real<N> &n, int shift)
  {
    high_precision_real<N> result;
    result.negative() = n.negative();
    if (shift > 0)
      detail::raw_shiftleft(n, result, shift);
    else if (shift < 0)
      detail::raw_shiftright(n, result, -shift);
    else
      result = n;
    return result;
  }

  template <int N>
  high_precision_real<N> operator>>(const high_precision_real<N> &n, int shift)
  {
    return n << -shift;
  }

  template <int N>
  bool valid_precision(const high_precision_real<N> &n)
  {
    for (int i = 0; i < n.size() - 2; ++i)
    {
      if (n[i])
        return true;
    }
    // Ensure we have at least 64+32 = 96 bits
    // We must have something in the top 32-bits
    return n[n.size() - 2] & 0xffffffff00000000ull;
  }

  template <int N>
  bool valid_precision_for_inverse(const high_precision_real<N> &n)
  {
    for (int i = 0; i < n.size() - 2; ++i)
    {
      if (n[i])
        return true;
    }
    // Ensure we have at least 64+48 = 112 bits
    // We must have something in the top 48-bits
    return n[n.size() - 2] & 0xffffffffff000000ull;
  }

  template <int N, int M>
  struct number_cast_t<
      high_precision_real<N>, high_precision_real<M>>
  {
    static high_precision_real<N> cast(const high_precision_real<M> &src)
    {
      return src;
    }
  };

  template <int N, native_floating_point From>
  struct number_cast_t<
      high_precision_real<N>, From>
  {
    static high_precision_real<N> cast(From src)
    {
      return {src};
    }
  };

  template <int N, non_native_floating_point From>
  struct number_cast_t<
      high_precision_real<N>, From>
  {
    static high_precision_real<N> cast(From src)
    {
      auto p = mantissa_exponent(src);
      return high_precision_real<N>{p.first}<<p.second;
    }
  };


  template <native_floating_point To, int N>
  struct number_cast_t<To, 
      high_precision_real<N>>
  {
    static To cast(const high_precision_real<N> &src)
    {
      return to_double(src);
    }
  };

  template<int N>
  double to_double(const high_precision_real<N> &src) { return src.to_double(); }

  template <int N>
  std::pair<double, int> mantissa_exponent(const high_precision_real<N> &x)
  {
    if (x[0] == 0)
    {
      int e = count_fractional_zeros(x);
      return {(double)(x << e).to_double(), -e};
    }
    else
    {
      return {x.to_double(), 0};
    }
  }

  template <int N>
  std::uint64_t *begin(high_precision_real<N> &x)
  {
    return x.begin();
  }

  template <int N>
  const std::uint64_t *begin(const high_precision_real<N> &x)
  {
    return x.begin();
  }

  template <int N>
  std::uint64_t *end(high_precision_real<N> &x)
  {
    return x.end();
  }

  template <int N>
  const std::uint64_t *end(const high_precision_real<N> &x)
  {
    return x.end();
  }

  template <int N>
  std::size_t size(const high_precision_real<N> &x)
  {
    return x.size();
  }

  template <int N>
  std::size_t int_size(const high_precision_real<N> &x)
  {
    return 1;
  }

  template <int N>
  double log(const high_precision_real<N> &n)
  {
    auto p = mantissa_exponent(n);
    return std::log(p.first) + std::log(2) * p.second;
  }

  template <int N>
  struct number_traits<high_precision_real<N>>
  {
    static constexpr bool is_number = true;
    static constexpr bool is_native = false;
    static constexpr bool is_floating_point = false;
    static constexpr bool is_real = true;
    using value_type = typename high_precision_real<N>::part_type;
  };
}

static_assert(numbers::number<numbers::high_precision_real<512>>);
static_assert(numbers::real<numbers::high_precision_real<512>>);
static_assert(numbers::fixed_point<numbers::high_precision_real<512>>);
