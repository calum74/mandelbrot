#pragma once
#include "rendering_sequence.hpp"
#include "view_pixmap.hpp"

namespace fractals {
class fractal_calculation;

class calculation_pixmap : public async_rendering_sequence {
public:
  calculation_pixmap(view_pixmap &pm, int stride, fractal_calculation &calculation);

  std::atomic<double> min_depth, max_depth;
  std::atomic<std::uint64_t> points_calculated;
  virtual void layer_complete2(int stride) = 0;

protected:
  fractal_calculation &calculation;
  view_pixmap &pixels;

  void calculate_point(int x, int y, int w) override;
  void interpolate_region(int cx, int cy, int x0, int y0, int x1, int y1);
  bool maybe_fill_region(int x0, int y0, int x1, int y1);
  void interpolate_region_smooth(int x0, int y0, int x1, int y1);

  void layer_complete(int stride, std::atomic<bool> &stop) override;
};

} // namespace fractals
