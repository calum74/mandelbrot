#include "calculation_pixmap.hpp"
#include "fractal_calculation.hpp"
#include <cassert>

fractals::calculation_pixmap::calculation_pixmap(
    pixmap_type &pm, int stride, fractal_calculation &calculation)
    : pixels(pm), async_rendering_sequence(pm.width(), pm.height(), stride),
      calculation(calculation) {

  min_depth = 0;
  max_depth = 0;
  points_calculated = 0;
}

void fractals::calculation_pixmap::calculate_point(int x, int y, int stride) {

  double depth;
  if (pixels(x, y).error == 0) {
    depth = pixels(x, y).value;
  } else {
    ++points_calculated;
    depth = calculation.calculate(x, y);
    pixels(x, y) = {depth, 0};
  }

  if (depth > 0) {
    // Technically this is a race condition
    // but we want to capture this here (and not in layer_complete())
    // because we need this in case we abort computation before the first
    // layer is complete.
    if (depth < min_depth || min_depth == 0)
      min_depth = depth;
    if (depth > max_depth)
      max_depth = depth;
  }

  if (stride > 1) {
    // Interpolate the region
    if (x > 0 && y > 0) {
      maybe_fill_region(x - stride, y - stride, x, y);
      interpolate_region_smooth(x - stride, y - stride, x, y);
    }

#if 0 
    auto d = stride / 2;
    int x0 = x - d;
    int x1 = x + d;
    int y0 = y - d;
    int y1 = y + d;
    if (x0 < 0)
      x0 = 0;
    if (x1 >= pixels.width())
      x1 = pixels.width() - 1;
    if (y0 < 0)
      y0 = 0;
    if (y1 >= pixels.height())
      y1 = pixels.height() - 1;
    interpolate_region(x, y, x0, y0, x1, y1);
#endif
  }
}

bool fractals::calculation_pixmap::maybe_fill_region(int x0, int y0, int x1,
                                                     int y1) {
  if (x1 - x0 > 2)
    return false;
  auto c00 = pixels(x0, y0);
  auto c10 = pixels(x1, y0);
  auto c01 = pixels(x0, y1);
  auto c11 = pixels(x1, y1);

  // If all 4 corners have the same colour, claim that the filled in colour is
  // accurate and does not need to be recalculated 1 means more speed 0 means
  // more accuracy
  if (c00.value == c10.value && c00.value == c11.value &&
      c00.value == c01.value) {
    for (int j = y0; j <= y1; ++j)
      for (int i = x0; i <= x1; ++i)
        pixels(i, j) = c00;
    return true;
  }
  return false;
}

void fractals::calculation_pixmap::interpolate_region_smooth(int x0, int y0,
                                                             int x1, int y1) {
  auto c00 = pixels(x0, y0);
  auto c10 = pixels(x1, y0);
  auto c01 = pixels(x0, y1);
  auto c11 = pixels(x1, y1);

  for (int j = y0; j <= y1; ++j) {
    for (int i = x0; i <= x1; ++i) {
      auto &p = pixels(i, j);
      int new_error = std::min(i - x0, x1 - i) + std::min(j - y0, y1 - j);

      if (new_error > 0 && new_error < p.error) {
        double px = (i - x0) / (double)(x1 - x0);
        double py = (j - y0) / (double)(y1 - y0);
        double new_value = c00.value * (1 - px) * (1 - py) +
                           c10.value * px * (1 - py) +
                           c01.value * (1 - px) * py + c11.value * px * py;
        p = {new_value, new_error};
      }
    }
  }
}

void fractals::calculation_pixmap::interpolate_region(int cx, int cy, int x0,
                                                      int y0, int x1, int y1) {
  auto &c = pixels(cx, cy);
  assert(x0 >= 0);
  assert(x1 < pixels.width());
  assert(y0 >= 0);
  assert(y0 < pixels.height());

  // Solid colour
  for (int j = y0; j <= y1; ++j) {
    int dy = j - cy;
    int ey = dy < 0 ? -dy : dy;
    for (int i = x0; i <= x1; ++i) {
      int dx = i - cx;
      int ex = dx < 0 ? -dx : dx;
      int error = ex + ey; // Manhattan distance
      auto &p = pixels(i, j);
      // This gives the artistic effects on the zoom
      if (error < p.error) {
        p = c;
        p.error = error;
      }
    }
  }
}

double fractals::fractal_calculation::average_iterations() const { return 0; }

double fractals::fractal_calculation::average_skipped() const { return 0; }

void fractals::fractal_calculation::initialize(const view_coords &c, int x,
                                               int y, std::atomic<bool> &stop) {
}
