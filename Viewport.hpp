#pragma once
#include "RGB.hpp"

namespace fractals {
/*
A region to render. It is essentially a buffer of pixels
*/
struct Viewport {
  int width, height;

  // The buffer to receive the pixels, stored in row-order.
  // data must be width*height elements large.
  RGB *data;

  RGB &operator()(int x, int y) { return data[x + y * width]; }
  RGB operator()(int x, int y) const { return data[x + y * width]; }

  // Callback when data has changed in the given pixel-region
  virtual void region_updated(int x, int y, int w, int h);

  // Callback when viewport fully calculated.
  // May not get called if cancelled.
  virtual void finished(double width, int min_depth, int max_depth,
                        double render_time);
};

} // namespace fractals