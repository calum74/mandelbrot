#pragma once
#include "pixmap.hpp"
#include "rendering_sequence.hpp"

namespace fractals {
class fractal_calculation;

class calculation_pixmap : public async_rendering_sequence {
public:
  using value_type = error_value<double>;
  using pixmap = pixmap<value_type>;

  calculation_pixmap(pixmap &pm, int stride, fractal_calculation &calculation);

  std::atomic<double> min_depth, max_depth;
  std::atomic<std::uint64_t> points_calculated;

protected:
  fractal_calculation &calculation;
  pixmap &pixels;

  void calculate_point(int x, int y, int w) override;
  void interpolate_region(int cx, int cy, int x0, int y0, int x1, int y1);
  bool maybe_fill_region(int x0, int y0, int x1, int y1);
};

} // namespace fractals
