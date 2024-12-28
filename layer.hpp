#pragma once

#include "view_coords.hpp"
#include <vector>

namespace fractals {

class PointwiseCalculation;

class layer {

public:
  layer(const view_coords &vc, int w, int h, const PointwiseCalculation &,
        std::atomic<bool> &stop, int threads);

  float &operator()(int x, int y) { return points[x + y * width]; }
  float operator()(int x, int y) const { return points[x + y * width]; }

private:
  view_coords coords;
  float calculation_time_seconds;
  int width, height;
  std::vector<float> points;
};
} // namespace fractals
