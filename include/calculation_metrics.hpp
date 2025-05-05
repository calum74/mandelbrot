#pragma once
#include "radius.hpp"

#include <cstddef>
#include <cstdint>

namespace fractals {
// For reporting on what happened during rendering a view
struct calculation_metrics {
  bool last_action_was_a_scroll;

  // These metrics are always reported
  numbers::radius radius;
  double min_depth;
  double max_depth;
  double p999, p9999;
  int p9999_x, p9999_y;
  std::uint64_t points_calculated;
  std::uint64_t non_black_points;
  double seconds_per_point;

  double discovered_depth; // The depth of the 99th percentile pixel

  bool fully_evaluated; // True if the calculation evaluated every point

  // These metrics are reported only on successful completion
  double average_iterations;
  double average_skipped_iterations;
  double render_time_seconds;
};
} // namespace fractals
