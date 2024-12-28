#pragma once
#include "high_precision_real.hpp"

namespace fractals {
struct view_coords {
  using value_type = high_precision_real<20>;
  value_type x, y, r;
  int max_iterations;

  view_coords scroll(int w, int h, int dx, int dy) const;

  view_coords zoom(double ratio, int w, int h, int cx, int cy) const;
};
} // namespace fractals
