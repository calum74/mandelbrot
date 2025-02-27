#pragma once
#include "calculation_metrics.hpp"
#include "fractal.hpp"
#include "pixmap.hpp"
#include "view_coords.hpp"
#include <atomic>
#include <future>

namespace fractals {
/*
This is currently a work in progress.

Manages the multithreaded calculation and animation
*/
class view {
public:
  view();
  ~view();

  pixmap<double> values;

  class listener {
  public:
    virtual void update(const calculation_metrics &) = 0;
    virtual void animation_timeout(const calculation_metrics &) = 0;
  };

  void set_size(int w, int h);
  void set_fractal(const fractal &new_fractal);
  void set_listener(listener *);
  void set_coords(const view_coords &new_coords);
  void set_threading(int threads); // 0 for hardware

  // Performs a smooth zoom in or a zoom out to the specified point.
  void animate_to(int x, int y, double ratio, double duration);

  // Performs an instant zoom in or out to the given point
  void zoom(int x, int y, double ratio);

  // Performs an instant scroll to the given point
  void scroll(int dx, int dy);

  void stop_current_animation_and_set_as_current();

private:
  std::shared_ptr<fractal_calculation> calculation;
  std::shared_ptr<fractal_calculation_factory> fractal;
  listener *current_listener;

  view_coords current_coords;  // What we're currently calculating
  calculation_metrics metrics;

  std::future<void> calculation_future;
  std::future<void> animation_future;
  std::atomic<bool> stop_calculation;
  std::atomic<bool> stop_animation;
  std::condition_variable animation_condition;

  pixmap<error_value<double>> previous_calculation_values, current_calculation_values;

  class my_calculation_pixmap;

  int calculation_threads;

  void start_calculating();
  void stop_calculating();
  void start_animating();
  void stop_animating();

  void animation_thread();

  void complete_layer(double min_depth, double max_depth, std::uint64_t points_calculated, int stride);
};
} // namespace fractals
