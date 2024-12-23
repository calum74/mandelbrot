#pragma once
#include "high_precision_real.hpp"

namespace fractals {
struct ViewCoords {
  using value_type = high_precision_real<20>;
  value_type x, y, r;
  int max_iterations;
};
} // namespace fractals
