/*
  Additional utilities and algorithms for dealing with complex numbers.
*/

#pragma once

#include "real_number.hpp"  // ???
#include "high_exponent_real.hpp"  // !! No

#include "number_traits.hpp"
#include "number_cast.hpp"


#include <complex>

namespace numbers
{
  template<typename T>
  T real_part(const std::complex<T> &c) { return c.real(); }

  template<typename T>
  T imag_part(const std::complex<T> &c) { return c.imag(); }
  
  template <real R>
  struct number_traits<std::complex<R>>
  {
    static constexpr bool is_number = true;
    static constexpr bool is_native = false;
    static constexpr bool is_floating_point = false;
    static constexpr bool is_real = false;
  };

  template<typename T>
  concept complex = // number<T> and not real<T> and
  requires(T t)
  {
      real_part(t);
      imag_part(t);
  };
}

namespace numbers
{

  template <typename T>
  concept Complex = requires(T v) {
    { v.real() } -> std::same_as<typename T::value_type>;
    { v.imag() } -> std::same_as<typename T::value_type>;
  };

  //template <typename T>
  //T real_part(const std::complex<T> &c) { return c.real(); }

  //template <typename T>
  //T imag_part(const std::complex<T> &c) { return c.imag(); }

  template <Complex C>
  C square(const C &c)
  {
    return {real_part(c) * real_part(c) - imag_part(c) * imag_part(c),
            typename C::value_type(2) * real_part(c) * imag_part(c)};
  }

  template <typename T>
  std::complex<T> operator*(int k, const std::complex<T> &c)
  {
    return {T(k) * real_part(c), T(k) * imag_part(c)};
  }

  template <typename T>
  std::complex<T> mul(const std::complex<T> &a, const std::complex<T> &b)
  {
    return {real_part(a) * real_part(b) - imag_part(a) * imag_part(b),
            imag_part(a) * real_part(b) + real_part(a) * imag_part(b)};
  }

  // inline double to_double(double d) { return d; }
  inline double to_double(float f) { return f; }

  // Specialise this because std::complex does incompatible things
  template <typename D, typename E, bool N1, bool N2>
  std::complex<numbers::high_exponent_real<D, E, false>>
  operator*(std::complex<numbers::high_exponent_real<D, E, N1>> a,
            std::complex<numbers::high_exponent_real<D, E, N2>> b)
  {
    return {a.real() * b.real() - a.imag() * b.imag(),
            a.real() * b.imag() + a.imag() * b.real()};
  }

  // Specialise this because std::complex does incompatible things
  template <typename D, typename E, bool N1, bool N2>
    requires(N1 != N2)
  std::complex<high_exponent_real<D, E, false>>
  operator/(const std::complex<high_exponent_real<D, E, N1>> &a,
            const std::complex<high_exponent_real<D, E, N2>> &b)
  {
    auto d = b.real() * b.real() + b.imag() * b.imag();
    return {(a.real() * b.real() + a.imag() * b.imag()) / d,
            (a.imag() * b.real() - a.real() * b.imag()) / d};
  }

  template <typename D, typename E, bool N>
  std::complex<high_exponent_real<D, E, false>>
  operator/(const std::complex<high_exponent_real<D, E, N>> &a,
            const std::complex<high_exponent_real<D, E, N>> &b)
  {
    auto d = b.real() * b.real() + b.imag() * b.imag();
    return {(a.real() * b.real() + a.imag() * b.imag()) / d,
            (a.imag() * b.real() - a.real() * b.imag()) / d};
  }

  template <typename D, typename E, bool N1, bool N2>
  std::complex<high_exponent_real<D, E, false>>
  operator+(std::complex<high_exponent_real<D, E, N1>> a,
            std::complex<high_exponent_real<D, E, N2>> b)
  {
    return {a.real() + b.real(), a.imag() + b.imag()};
  }

  template <typename D, typename E, bool N>
  std::complex<high_exponent_real<D, E, false>>
  operator+(std::complex<high_exponent_real<D, E, N>> a,
            std::complex<high_exponent_real<D, E, N>> b)
  {
    return {a.real() + b.real(), a.imag() + b.imag()};
  }

  template <typename D, typename E, bool N>
  std::complex<high_exponent_real<D, E, false>>
  operator*(std::complex<high_exponent_real<D, E, N>> a,
            std::complex<high_exponent_real<D, E, N>> b)
  {
    return {a.real() * b.real() - a.imag() * b.imag(),
            a.real() * b.imag() + a.imag() * b.real()};
  }

  // template <Complex C> C step(const C &z, const C &c) { return square(z) + c; }

  template <Complex C>
  auto norm(const C &c)
  {
    auto r = real_part(c);
    auto i = imag_part(c);
    return r * r + i * i;
  }

  template <int Order, Complex C, bool is_even = (Order % 2 == 0)>
  struct pow_impl;

  template <Complex C>
  struct pow_impl<2, C, true>
  {
    static C eval(const C &c) { return square(c); }
  };

  template <Complex C>
  struct pow_impl<1, C, false>
  {
    static C eval(const C &c) { return c; }
  };

  template <int Order, Complex C>
  struct pow_impl<Order, C, true>
  {
    static C eval(const C &c)
    {
      auto r = pow_impl<Order / 2, C>::eval(c);
      return square(r);
    }
  };

  template <Complex C>
  struct pow_impl<0, C, true>
  {
    static C eval(const C &c) { return C{1, 0}; }
  };

  template <Complex C>
  struct pow_impl<-1, C, false>
  {
    static C eval(const C &c) { return C{0, 0}; } // Not implemented
  };

  template <int Order, Complex C>
  struct pow_impl<Order, C, false>
  {
    static C eval(const C &c) { return mul(c, pow_impl<Order - 1, C>::eval(c)); }
  };

  template <int Order, Complex C>
  C pow(const C &c)
  {
    return pow_impl<Order, C>::eval(c);
  }

  template <int N, int M>
  struct choose_impl
  {
    static const int value =
        choose_impl<N - 1, M - 1>::value + choose_impl<N - 1, M>::value;
  };

  template <int M>
  struct choose_impl<0, M>
  {
    static const int value = 0;
  };

  template <>
  struct choose_impl<0, 1>
  {
    static const int value = 0;
  };

  template <int N>
  struct choose_impl<N, N>
  {
    static const int value = 1;
  };

  template <int N>
  struct choose_impl<N, 0>
  {
    static const int value = 1;
  };

  template <int N>
  struct choose_impl<N, 1>
  {
    static const int value = N;
  };

  template <int N, int M>
  constexpr int choose()
  {
    return choose_impl<N, M>::value;
  }

  inline bool isfinite(double d) { return std::isfinite(d); }

  template <typename T>
  std::complex<T> normalize(const std::complex<T> &c)
  {
    return {normalize(c.real()), normalize(c.imag())};
  }

  template <typename T>
  struct normalized<std::complex<T>>
  {
    using type = std::complex<typename normalized<T>::type>;
  };

  template <real T, typename From>
  struct number_cast_t<std::complex<T>, From>
  {
    static std::complex<T> cast(const From &x)
    {
      return {number_cast<T>(real_part(x)), number_cast<T>(imag_part(x))};
    }
  };

  // !! Rename to complex_generator

  // A complex number with the given precision
  template <int Digits, int MinExp, int MaxExp>
  using complex_number =
      std::complex<typename make_real<Digits, MinExp, MaxExp>::type>;
} // namespace fractals


static_assert(numbers::complex<std::complex<double>>);