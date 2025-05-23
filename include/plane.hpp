#pragma once
#include "view_coords.hpp"
#include "numbers/number_cast.hpp"

namespace fractals {

template <typename T, typename Tag = void> struct point {
  T x, y;
};

// A plane maps points in the display to points in the fractal.
template <typename C1, typename C2> class plane {
  public:
  C1 x0, y0;
  C2 w, h, dx, dy;
  int pw, ph;

  plane() = default;

  plane(const view_coords &c, int pw, int ph) : pw(pw), ph(ph) {
    if (pw > ph) {
      y0 = number_cast<C1>(c.y) - number_cast<C1>(c.r);
      h = number_cast<C2>(c.r) * C2(2);
      w = number_cast<C2>(c.r) * C2((2.0 * pw) / ph);
      x0 = number_cast<C1>(c.x) - numbers::number_cast<C1>(w * C2(0.5));
    } else {
      x0 = number_cast<C1>(c.x) - number_cast<C1>(c.r);
      w = number_cast<C2>(c.r) * C2(2);
      h = number_cast<C2>(c.r) * C2((2.0 * ph) / pw);
      y0 = number_cast<C1>(c.y) - numbers::number_cast<C1>(h * C2(0.5));
    }
    dx = w * C2(1.0 / pw);
    dy = h * C2(1.0 / ph);
  }

  C1 get_x(int x) const { return x0 + numbers::number_cast<C1>(x * dx); }
  C1 get_y(int y) const { return y0 + numbers::number_cast<C1>(y * dy); }

  point<C1> operator()(int x, int y) const { return {get_x(x), get_y(y)}; }

  // Returns -1 if the point is out of range
  int to_x(const C1 &c) const
  {
    int x = numbers::number_cast<int>(numbers::number_cast<C2>(c - x0)/dx);
    return x>=0 && x<pw ? x : -1;
  }

  int to_y(const C1 &c) const
  {
    int y = numbers::number_cast<int>(numbers::number_cast<C2>(c - y0)/dy);
    return y>=0 && y<ph ? y : -1;
  }

};
} // namespace fractals
