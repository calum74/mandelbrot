#pragma once
#include "high_precision_real.hpp"

namespace fractals {

// The coordinates of the current view, specified in high_precision numbers.
struct view_coords {
  using value_type = high_precision_real<20>;
  value_type x, y, r;
  int max_iterations;

  // Scrolls the coordinates by a given number of pixels, returning the
  // resulting coordinates.
  view_coords scroll(int w, int h, int dx, int dy) const;

  // Zooms the coordinates in or out by the given ratio, returning the resulting
  // coordinates. ratio is the new size/old size, so ratio < 1 means to zoom in.
  // The zoom is centered around the point (cx,cy).
  view_coords zoom(double ratio, int w, int h, int cx, int cy) const;

  double point_size(int w, int h) const;
  value_type top(int w, int h) const;
  value_type left(int w, int h) const;
};

std::ostream &operator<<(std::ostream &os, const view_coords &coords);

} // namespace fractals
