#pragma once

/*
    Biassed floating-point numbers.
*/

#include "real_number.hpp"
#include <cmath>
#include <iostream>

namespace fractals {

struct raw_value {};

/*
    A biassed real is a floating point number whose real exponent is biassed by
   a fixed quantity `Bias`. This is typically used for numbers that have a very
   large exponent that cannot be represented by standard hardware numbers. But
   we can use templates to exploit
*/
template <typename Real, int Bias> class biassed_real {
public:
  // This may lose precision
  // explicit biassed_real(Real r) : value(r * std::exp2(-Bias)) {}

  biassed_real(Real r, raw_value) : value(r) {}

  template <int Bias2>
  explicit biassed_real(
      biassed_real<Real, Bias2> other) /* require (no loss of precision) */
      : value(other.value * std::exp2(Bias2 - Bias)) {}

  // Will likely lose precision and just not work
  double to_double() const { return value * std::exp2(Bias); }

  Real value;
};

template <typename Real, int Bias>
biassed_real<Real, Bias> operator+(biassed_real<Real, Bias> a,
                                   biassed_real<Real, Bias> b) {
  return {a.value + b.value, raw_value()};
}

template <typename Real, int Bias1, int Bias2>
biassed_real<Real, Bias1> operator+(biassed_real<Real, Bias1> a,
                                    biassed_real<Real, Bias2> b)
  requires(
      Bias2 + std::numeric_limits<Real>::max_exponent <=
      Bias1 +
          std::numeric_limits<Real>::min_exponent) // The numbers cannot overlap
{
  return a;
}

template <typename Real, int Bias1, int Bias2>
biassed_real<Real, Bias1> operator+(biassed_real<Real, Bias1> a,
                                    biassed_real<Real, Bias2> b)
  requires(Bias1 > Bias2 + std::numeric_limits<Real>::digits)
{
  return b;
}

template <typename Real, int Bias>
biassed_real<Real, Bias> operator-(biassed_real<Real, Bias> a,
                                   biassed_real<Real, Bias> b) {
  return {a.value - b.value, raw_value()};
}

template <typename Real, int Bias1, int Bias2>
biassed_real<Real, Bias1 + Bias2> operator*(biassed_real<Real, Bias1> a,
                                            biassed_real<Real, Bias2> b) {
  return {a.value + b.value, raw_value()};
}

template <typename Real, int Bias>
std::ostream &operator<<(std::ostream &os, const biassed_real<Real, Bias> &r) {
  return os << r.value << "*2^" << Bias;
}

template <typename Real, int Bias1, int Bias2>
std::complex<biassed_real<Real, Bias1 + Bias2>>
operator*(std::complex<Real, Bias1>, std::complex<Real, Bias2>);

} // namespace fractals
