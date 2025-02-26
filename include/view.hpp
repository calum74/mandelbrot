#pragma once
#include "pixmap.hpp"
#include "calculation_metrics.hpp"
#include "fractal.hpp"
#include "view_coords.hpp"
#include <atomic>
#include <future>

namespace fractals {
// Manages the multithreaded animation and rendering
class view {
public:
  pixmap<error_value<double>> values;

  class listener {
  public:
    virtual void redraw() = 0;
    virtual void update_status(const calculation_metrics &) = 0;
  };

  void set_size(int w, int h);
  void set_fractal(const fractal &new_fractal);
  void set_listener(listener &);
  void set_coords(const view_coords &new_coords);

private:
  std::shared_ptr<fractal_calculation> calculation;
  std::shared_ptr<fractal_calculation_factory> fractal;

  view_coords current_coords;
  calculation_metrics metrics;

  std::future<void> calculation_thread;
  std::future<void> animation_thread;
  std::atomic<bool> cancel_calculation;
  std::atomic<bool> cancel_animation;
};
} // namespace fractals
