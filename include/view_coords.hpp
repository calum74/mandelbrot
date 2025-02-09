#pragma once
#include "high_precision_real.hpp"

namespace fractals {
struct view_parameters;

struct mapped_point {
  // The coordinates of a point (x,y) within a viewport, together with
  // the log of the distance log(r1/r2).
  double x, y, log_distance;
};

// The coordinates of the current view, specified in high_precision numbers.
struct view_coords {
  using value_type = real_number<4096, 0, 0>; // 4096 with no exponent
  value_type x, y, r;
  int max_iterations;

  view_coords() = default;
  view_coords(const value_type &x, const value_type &y, const value_type &r,
              int max_iterations);
  view_coords(const view_parameters &);

  void write(view_parameters&) const;

  // Gets the natural log of r
  double ln_r() const;

  // Scrolls the coordinates by a given number of pixels, returning the
  // resulting coordinates.
  view_coords scroll(int w, int h, int dx, int dy) const;

  // Zooms the coordinates in or out by the given ratio, returning the resulting
  // coordinates. ratio is the new size/old size, so ratio < 1 means to zoom in.
  // The zoom is centered around the point (cx,cy).
  view_coords zoom(double ratio, int w, int h, int cx, int cy) const;

  view_coords zoom(double ratio, int w, int h, int cx, int cy,
                   const value_type &CX, const value_type &CY) const;

  // Zoom in on the center
  view_coords zoom(double ratio) const;

  double point_size(int w, int h) const;
  value_type point_size_full(int w, int h) const;
  value_type top(int w, int h) const;
  value_type left(int w, int h) const;

  std::pair<value_type, value_type> map_point(int w, int h, int cx,
                                              int cy) const;

  mapped_point map_point(int w, int h, const view_coords &p) const;

  // Gets the desired precision
  int get_precision(int dps = 4) const;
};

std::ostream &operator<<(std::ostream &os, const view_coords &coords);
std::istream &operator>>(std::istream &is, view_coords &coords);

// Helper to output the radius from its log in engineering format
void log_radius(std::ostream &os, double log_base_e);

} // namespace fractals
