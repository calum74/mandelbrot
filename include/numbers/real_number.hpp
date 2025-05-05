#pragma once
#include <limits>
#include "high_precision_real.hpp"
#include "high_exponent_real.hpp"

/*
  Selects the datatype that is able to represent a given type requirement.

  Defines real numbers based on desired characteristics.
  Instead of saying "I want a double", we can say "I want a number with
  precision X".

  This is a template library that finds the most suitable representation of a
  number, and can make use of any platform hardware which may exist on your
  platform. E.g. long double is not available on ARM64.
*/

namespace numbers
{
  template <int Digits, int MinExp, int MaxExp>
  struct make_real;

  template <int Digits, int MinExp, int MaxExp>
    requires(Digits <= std::numeric_limits<float>::digits &&
             MinExp >= std::numeric_limits<float>::min_exponent &&
             MaxExp <= std::numeric_limits<float>::max_exponent)
  struct make_real<Digits, MinExp, MaxExp>
  {
    using type = float;
  };

  template <int Digits, int MinExp, int MaxExp>
    requires((Digits > std::numeric_limits<float>::digits || MinExp<std::numeric_limits<float>::min_exponent || MaxExp> std::numeric_limits<float>::max_exponent) &&
             Digits <= std::numeric_limits<double>::digits &&
             MinExp >= std::numeric_limits<double>::min_exponent and
             MaxExp <= std::numeric_limits<double>::max_exponent)
  struct make_real<Digits, MinExp, MaxExp>
  {
    using type = double;
  };

  template <int Digits, int MinExp, int MaxExp>
    requires(Digits > std::numeric_limits<double>::digits &&
             Digits <= std::numeric_limits<long double>::digits &&
             MinExp >= std::numeric_limits<long double>::min_exponent and
             MaxExp <= std::numeric_limits<long double>::max_exponent)
  struct make_real<Digits, MinExp, MaxExp>
  {
    using type = long double;
  };

  inline bool valid_precision_for_inverse(double d) { return d >= 1e5 * std::numeric_limits<double>::epsilon(); }

  inline bool valid_precision(double d) { return d >= 1e7 * std::numeric_limits<double>::epsilon(); }

  inline double inverse(double d) { return 1 / d; }

  inline double normalize(double d) { return d; }
  inline float normalize(float f) { return f; }

  template <int Digits, int MinExp, int MaxExp>
  using real_number = typename make_real<Digits, MinExp, MaxExp>::type;

  template <typename T>
  struct normalized
  {
    using type = T;
  };

  template <int Digits>
    requires(Digits > std::numeric_limits<long double>::digits)
  struct make_real<Digits, 0, 0>
  {
    using type = high_precision_real<Digits>;
  };

  template <typename Double, typename Exponent>
  struct normalized<high_exponent_real<Double, Exponent, false>>
  {
    using type = high_exponent_real<Double, Exponent, true>;
  };

  template <int Digits, int MinExp, int MaxExp>
    requires(Digits > std::numeric_limits<double>::digits && Digits <= std::numeric_limits<long double>::digits &&
             (MinExp<std::numeric_limits<long double>::min_exponent || MaxExp>
                  std::numeric_limits<long double>::max_exponent))
  struct make_real<Digits, MinExp, MaxExp>
  {
    using type = high_exponent_real<long double, int>;
  };

  template <int Digits, int MinExp, int MaxExp>
    requires(Digits <= std::numeric_limits<double>::digits &&
             (MinExp<std::numeric_limits<double>::min_exponent || MaxExp>
                  std::numeric_limits<double>::max_exponent))
  struct make_real<Digits, MinExp, MaxExp>
  {
    using type = high_exponent_real<double, int>;
  };

}