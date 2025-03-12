#pragma once
#include "calculation_metrics.hpp"
#include "fractal.hpp"
#include "pixmap.hpp"
#include "radius.hpp"
#include "view_coords.hpp"
#include "view_pixmap.hpp"
#include <atomic>
#include <future>
#include <chrono>

namespace fractals {
class view_listener;

/*
This is currently a work in progress.

Manages the multithreaded calculation and animation
*/
class view {
public:
  view();
  ~view();

  view_pixmap values;

  void set_size(int w, int h);
  void set_fractal(const fractal &new_fractal, bool init_coords, bool recalculate);
  void set_listener(view_listener *);
  void set_coords(const view_coords &new_coords, bool recalculate);
  void set_threading(int threads); // 0 for hardware

  // Performs a smooth zoom x2 to the specified point.
  void animate_to(int x, int y, std::chrono::duration<double> duration,
                  bool wait_for_completion);
  void animate_to_center(std::chrono::duration<double> duration,
                         bool wait_for_completion, double ratio);

  // Performs an instant zoom in or out to the given point
  // Cancels any animations in progress
  void zoom(int x, int y, double ratio);

  // Performs an instant scroll to the given point
  // Cancels any animations in progress
  void scroll(int dx, int dy);

  void stop_current_animation_and_set_as_current();
  void start_calculating();

  const view_coords &get_coords() const;
  const calculation_metrics &get_metrics() const;

  void increase_iterations();
  void decrease_iterations();
  void set_max_iterations(int n);
  bool get_auto_zoom(int &x, int &y);
  void update_iterations(const calculation_metrics &);
  std::string get_fractal_name() const;
  std::string get_fractal_family() const;
  view_coords initial_coords() const;

  int width() const { return values.width(); }
  int height() const { return values.height(); }
  bool is_animating() const;

private:
  std::shared_ptr<fractal_calculation> calculation;
  std::shared_ptr<fractal_calculation_factory> fractal;
  view_listener *listener;

  view_coords calculation_coords; // What we're currently calculating

  std::future<void> calculation_future;
  std::future<void> animation_future;
  std::atomic<bool> stop_calculation;
  std::atomic<bool> stop_animation;

  std::mutex mutex;
  std::condition_variable animation_condition;
  std::atomic<bool> calculation_completed; // Protected by mutex
  calculation_metrics metrics;             // Protected by mutex

  // True if we are calculating in the background.
  // This means that we'll display a zoom of the background image,
  // not the currently calculated image.
  bool animating;  // Protected by mutex

  // True in quality rendering mode,
  // where we only display the new image once it has completely finished
  // calculating.
  bool wait_for_calculation_to_complete;
  bool paused_waiting_for_calculation;

  // Which point we are zooming into
  int zoom_x, zoom_y;

  // How far we have zoomed so far
  // 1.0 = just started
  // 0.5 = finished zooming
  double rendered_zoom_ratio;

  std::chrono::time_point<std::chrono::system_clock> animation_start;
  std::chrono::duration<double> animation_duration;

  view_pixmap previous_calculation_values, current_calculation_values;

  class my_calculation_pixmap;

  int calculation_threads;

  void stop_calculating();
  void start_animating();
  void stop_animating();

  void animation_thread();
  void complete_layer(double min_depth, double max_depth,
                      std::uint64_t points_calculated, int stride);
  bool valid() const;

  void freeze_current_view();
  void animate_to(int x, int y, std::chrono::duration<double> duration,
    bool wait_for_completion, bool lock_center, double ratio);

};

void measure_depths(const view_pixmap & values, calculation_metrics &metrics);

} // namespace fractals
