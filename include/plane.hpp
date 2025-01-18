#pragma once
#include "view_coords.hpp"

namespace fractals {

template <typename T, typename Tag = void> struct point {
  T x, y;
};

// A plane maps points in the display to points in the fractal.
template <typename C1, typename C2 = C1> struct plane {
  C1 x0, y0;
  C2 w, h, dx, dy;

  plane() = default;

  plane(const view_coords &c, int pw, int ph) {
    if (pw > ph) {
      y0 = convert<C1>(c.y) - convert<C1>(c.r);
      h = convert<C2>(c.r) * C2(2);
      w = convert<C2>(c.r) * C2((2.0 * pw) / ph);
      x0 = convert<C1>(c.x) - convert<C1>(w * C2(0.5));
    } else {
      x0 = convert<C1>(c.x) - convert<C1>(c.r);
      w = convert<C2>(c.r) * C2(2);
      h = convert<C2>(c.r) * C2((2.0 * ph) / pw);
      y0 = convert<C1>(c.y) - convert<C1>(h * C2(0.5));
    }
    dx = w * C2(1.0 / pw);
    dy = h * C2(1.0 / ph);
  }

  C1 get_x(int x) const { return x0 + x * dx; }
  C1 get_y(int y) const { return y0 + y * dy; }

  point<C1> operator()(int x, int y) const { return {get_x(x), get_y(y)}; }
};
} // namespace fractals
