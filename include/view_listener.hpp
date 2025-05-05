#pragma once
#include "radius.hpp"

namespace fractals {
class calculation_metrics;

// Called back on any thread, but will not overlap
class view_listener {
public:
  virtual ~view_listener() = default;
  // A new calculation has started
  virtual void calculation_started(numbers::radius r, int max_iterations) = 0;

  // The contents of `values` has changed, so it can be redrawn
  // Triggered by animation or by calculation.
  virtual void values_changed() = 0;

  // A calculation has completed or was interrupted.
  // Called once per calculation.
  virtual void calculation_finished(const calculation_metrics &) = 0;

  // Called when it's time to schedule the next animation step.
  // The calculation may or may not have already finished.
  // If in `wait` mode, will only be called after `calculation_finished`.
  virtual void animation_finished(const calculation_metrics &) = 0;
};

} // namespace fractals