#include "view_pixmap.hpp"

#include <numbers>

constexpr fractals::error_value<double> missing_value{
    std::numeric_limits<double>::quiet_NaN(), 127};

namespace fractals {
template <typename ErrorFn>
void interpolate_values(const view_pixmap &src, view_pixmap &dest, double dx,
                        double dy, double r, ErrorFn fn) {

  for (int j = 0; j < dest.height(); ++j)
    for (int i = 0; i < dest.width(); ++i) {
      double rx = r * i + dx, ry = r * j + dy;
      int i2 = rx;
      int j2 = ry;
      auto &to_pixel = dest(i, j);
      if (i2 >= 0 && i2 < dest.width() && j2 >= 0 && j2 < dest.height()) {
        rx -= i2;
        ry -= j2;
        auto &from_pixel_00 = src(i2, j2);
        auto &from_pixel_10 = src(i2 + 1 < dest.width() ? i2 + 1 : i2, j2);
        auto &from_pixel_01 = src(i2, j2 + 1 < dest.height() ? j2 + 1 : j2);
        auto &from_pixel_11 = src(i2 + 1 < dest.width() ? i2 + 1 : i2,
                                  j2 + 1 < dest.height() ? j2 + 1 : j2);

        auto to_value = from_pixel_00.value * (1 - rx) * (1 - ry) +
                        from_pixel_10.value * rx * (1 - ry) +
                        from_pixel_01.value * (1 - rx) * ry +
                        from_pixel_11.value * rx * ry;
        auto to_error =
            std::max(std::max(from_pixel_00.error, from_pixel_10.error),
                     std::max(from_pixel_01.error, from_pixel_11.error));
        to_pixel = {to_value, fn(to_error)};
      } else {
        to_pixel = missing_value;
      }
    }
}
} // namespace fractals

void fractals::interpolate_values(const view_pixmap &src, view_pixmap &dest,
                                  double dx, double dy, double r) {
  interpolate_values(src, dest, dx, dy, r, [&](int e) { return e; });
}

void fractals::map_values(const view_pixmap &src, view_pixmap &dest, double dx,
                          double dy, double r) {
  bool zoom_eq = r == 1.0;
  bool zoom_out = r > 1.0;

  interpolate_values(src, dest, dx, dy, r, [&](int e) {
    return zoom_eq ? e : zoom_out ? 20 : e > 20 ? e : e + 1;
  });
}

void fractals::invalidate_values(view_pixmap &pm) {
  for (auto &p : pm)
    p.error = 127;
}
