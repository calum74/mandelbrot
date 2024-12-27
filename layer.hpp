#pragma once

#include "view_coords.hpp"
#include <vector>

namespace fractals {

// A layer is a rendered
class layer {

public:
  layer(int w, int h);
  float &operator()(int x, int y) { return points[x + y * width]; }
  float operator()(int x, int y) const { return points[x + y * width]; }

private:
  view_coords view;
  float calculation_time_seconds;
  int width, height;
  std::vector<float> points;
};
} // namespace fractals
