#include "view_coords.hpp"

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
  double pw = w;
  double ph = h;

  auto pixel_width = w > h ? r * (2.0 / ph) : r * (2.0 / pw);

  auto point_size = w > h ? r * (2.0 / ph) : r * (2.0 / pw);

  view_coords::value_type ratio2{ratio};

  auto CX = x + pixel_width * (cx - pw / 2);
  auto CY = y + pixel_width * (cy - ph / 2);

  view_coords new_coords;
  new_coords.max_iterations = max_iterations;
  new_coords.x = CX - (CX - x) * ratio2;
  new_coords.y = CY - (CY - y) * ratio2;
  new_coords.r = r * ratio2;

  return new_coords;
}