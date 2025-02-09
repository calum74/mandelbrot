#include "view_coords.hpp"
#include "high_exponent_real.hpp"
#include "view_parameters.hpp"
#include <iomanip>
#include <numbers>

fractals::view_coords::view_coords(const value_type &x, const value_type &y,
                                   const value_type &r, int max_iterations)
    : x(x), y(y), r(r), max_iterations(max_iterations) {}

fractals::view_coords::view_coords(const view_parameters &vp) {
  std::stringstream(vp.x) >> x;
  std::stringstream(vp.y) >> y;
  std::stringstream(vp.r) >> r;
  max_iterations = vp.max_iterations;
}

void fractals::view_coords::write(view_parameters&vp) const {
  auto precision = get_precision();
  std::stringstream sx;
  sx << std::setprecision(precision) << x;
  vp.x = sx.str();

  std::stringstream sy;
  sy << std::setprecision(precision) << y;
  vp.y = sy.str();

  std::stringstream sr;
  sr << std::setprecision(4);
  fractals::log_radius(sr, ln_r());
  vp.r = sr.str();

  vp.max_iterations = max_iterations;
}

fractals::view_coords fractals::view_coords::scroll(int w, int h, int dx,
                                                    int dy) const {
  auto result = *this;
  if (w > h) {
    result.y += r * (2.0 * dy / h);
    result.x += r * (2.0 * dx / h);
  } else {
    result.y += r * (2.0 * dy / w);
    result.x += r * (2.0 * dx / w);
  }
  return result;
}

fractals::view_coords fractals::view_coords::zoom(double ratio, int w, int h,
                                                  int cx, int cy) const {
  auto C = map_point(w, h, cx, cy);
  return zoom(ratio, w, h, cx, cy, C.first, C.second);
}

fractals::view_coords fractals::view_coords::zoom(double ratio) const {
  auto c = *this;
  c.r = c.r * ratio;
  return c;
}

std::pair<fractals::view_coords::value_type, fractals::view_coords::value_type>
fractals::view_coords::map_point(int w, int h, int cx, int cy) const {
  double pw = w;
  double ph = h;

  auto pixel_width = w > h ? r * (2.0 / ph) : r * (2.0 / pw);

  auto CX = x + pixel_width * (cx - pw / 2);
  auto CY = y + pixel_width * (cy - ph / 2);
  return {CX, CY};
}

fractals::view_coords fractals::view_coords::zoom(double ratio, int w, int h,
                                                  int cx, int cy,
                                                  const value_type &CX,
                                                  const value_type &CY) const {

  view_coords::value_type ratio2{ratio};
  view_coords new_coords;
  new_coords.max_iterations = max_iterations;
  new_coords.x = CX - (CX - x) * ratio2;
  new_coords.y = CY - (CY - y) * ratio2;
  new_coords.r = r * ratio2;

  return new_coords;
}

double fractals::view_coords::point_size(int w, int h) const {
  double pw = w;
  double ph = h;
  double pr = r.to_double();

  return pr * (w > h ? (2.0 / ph) : (2.0 / pw));
}

fractals::view_coords::value_type
fractals::view_coords::point_size_full(int w, int h) const {
  double pw = w;
  double ph = h;
  return r * (w > h ? (2.0 / ph) : (2.0 / pw));
}

fractals::view_coords::value_type fractals::view_coords::top(int w,
                                                             int h) const {

  return w > h ? y - r : y - r * (double(h) / double(w));
}

fractals::view_coords::value_type fractals::view_coords::left(int w,
                                                              int h) const {

  return w > h ? x - r * (double(w) / double(h)) : x - r;
}

std::ostream &fractals::operator<<(std::ostream &os,
                                   const view_coords &coords) {
  int width = coords.get_precision();

  return os << '(' << std::setprecision(width) << coords.x << ',' << coords.y
            << ',' << coords.r << ',' << coords.max_iterations << ')';
}

std::istream &fractals::operator>>(std::istream &is, view_coords &coords) {
  // !! This has no error detection and recovery whatsoever
  char ch1, ch2, ch3, ch4, ch5;
  return is >> ch1 >> coords.x >> ch2 >> coords.y >> ch3 >> coords.r >> ch4 >>
         coords.max_iterations >> ch5;
}

int fractals::view_coords::get_precision(int d) const {
  int zeros = fractals::count_fractional_zeros(r);
  return d + zeros * 0.30103;
}

void fractals::log_radius(std::ostream &os, double log_base_e) {
  // Renders a number ln(x) in engineering form
  auto log_base_10 = log_base_e * std::numbers::log10e;

  double int_part, frac_part = std::pow(10, std::modf(log_base_10, &int_part));
  while (frac_part < 1) {
    frac_part *= 10;
    int_part--;
  }
  os << frac_part << "e" << (int)int_part;
}

double fractals::view_coords::ln_r() const {
  return fractals::log(
      fractals::convert<fractals::high_exponent_real<double>>(r));
}

fractals::mapped_point
fractals::view_coords::map_point(int w, int h, const view_coords &p) const {

  // !! Use exponents here
  auto dx = (p.x - x).to_double();
  auto dy = (p.y - y).to_double();
  auto dr = r.to_double();
  auto point_size = (w > h ? double(h) / (2.0 * dr) : double(w) / (2.0 * dr));

  return {w / 2 + point_size * dx, h / 2 + point_size * dy, ln_r() - p.ln_r()};
}
