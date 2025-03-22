#include "radius.hpp"

#include <numbers>

fractals::radius::radius(double v, from_ln) : ln_r_value(v) {}

fractals::radius::radius() : ln_r_value(0) {}

fractals::radius::radius(double v) : ln_r_value(std::log(v)) {}

double fractals::radius::ln_r() const { return ln_r_value; }

double fractals::radius::to_double() const { return std::exp(ln_r_value); }

fractals::radius fractals::operator/(fractals::radius a, fractals::radius b) {
  return {a.ln_r() - b.ln_r(), radius::from_ln{}};
}

fractals::radius fractals::operator*(fractals::radius a, fractals::radius b) {
  return {a.ln_r() + b.ln_r(), radius::from_ln{}};
}

std::ostream &fractals::operator<<(std::ostream &os, radius r) {
  auto log_base_10 = r.ln_r() * std::numbers::log10e;

  double int_part, frac_part = std::pow(10, std::modf(log_base_10, &int_part));
  while (frac_part < 1) {
    frac_part *= 10;
    int_part--;
  }
  os << frac_part << "e" << (int)int_part;
  return os;
}

bool fractals::operator==(radius a, radius b) { return a.ln_r() == b.ln_r(); }
bool fractals::operator<(radius a, radius b) { return a.ln_r() < b.ln_r(); }
bool fractals::operator<=(radius a, radius b) { return a.ln_r() <= b.ln_r(); }
bool fractals::operator>(radius a, radius b) { return a.ln_r() > b.ln_r(); }
bool fractals::operator>=(radius a, radius b) { return a.ln_r() >= b.ln_r(); }
bool fractals::operator!=(radius a, radius b) { return a.ln_r() != b.ln_r(); }
