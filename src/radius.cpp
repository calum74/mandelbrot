#include "numbers/radius.hpp"

#include <numbers>
#include <cmath>

numbers::radius::radius(double v, from_ln) : ln_r_value(v) {}

numbers::radius::radius() : ln_r_value(0) {}

numbers::radius::radius(double v) : ln_r_value(std::log(v)) {}

double numbers::radius::ln_r() const { return ln_r_value; }

double numbers::radius::to_double() const { return std::exp(ln_r_value); }

numbers::radius numbers::operator/(numbers::radius a, numbers::radius b) {
  return {a.ln_r() - b.ln_r(), radius::from_ln{}};
}

numbers::radius numbers::operator*(numbers::radius a, numbers::radius b) {
  return {a.ln_r() + b.ln_r(), radius::from_ln{}};
}

std::ostream &numbers::operator<<(std::ostream &os, radius r) {
  auto log_base_10 = r.ln_r() * std::numbers::log10e;

  double int_part, frac_part = std::pow(10, std::modf(log_base_10, &int_part));
  while (frac_part < 1) {
    frac_part *= 10;
    int_part--;
  }
  os << frac_part << "e" << (int)int_part;
  return os;
}

bool numbers::operator==(radius a, radius b) { return a.ln_r() == b.ln_r(); }
bool numbers::operator<(radius a, radius b) { return a.ln_r() < b.ln_r(); }
bool numbers::operator<=(radius a, radius b) { return a.ln_r() <= b.ln_r(); }
bool numbers::operator>(radius a, radius b) { return a.ln_r() > b.ln_r(); }
bool numbers::operator>=(radius a, radius b) { return a.ln_r() >= b.ln_r(); }
bool numbers::operator!=(radius a, radius b) { return a.ln_r() != b.ln_r(); }
