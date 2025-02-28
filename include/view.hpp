#pragma once
#include "calculation_metrics.hpp"
#include "fractal.hpp"
#include "pixmap.hpp"
#include "view_coords.hpp"
#include <atomic>
#include <future>

namespace fractals {

using view_pixmap = pixmap<error_value<double>>;

/*
This is currently a work in progress.

Manages the multithreaded calculation and animation
*/
class view {
public:
  view();
  ~view();

  view_pixmap values;

  class listener {
  public:
    virtual void redraw() = 0;
    virtual void update(const calculation_metrics &) = 0;
    virtual void animation_complete(const calculation_metrics &) = 0;
  };

  void set_size(int w, int h);
  void set_fractal(const fractal &new_fractal);
  void set_listener(listener *);
  void set_coords(const view_coords &new_coords);
  void set_threading(int threads); // 0 for hardware

  // Performs a smooth zoom x2 to the specified point.
  void animate_to(int x, int y, std::chrono::duration<double> duration, bool wait_for_completion);
  void animate_to_center(std::chrono::duration<double> duration, bool wait_for_completion);

  // Performs an instant zoom in or out to the given point
  // Cancels any animations in progress
  void zoom(int x, int y, double ratio);

  // Performs an instant scroll to the given point
  // Cancels any animations in progress
  void scroll(int dx, int dy);

  void stop_current_animation_and_set_as_current();
  void start_calculating();

private:
  std::shared_ptr<fractal_calculation> calculation;
  std::shared_ptr<fractal_calculation_factory> fractal;
  listener *current_listener;

  view_coords current_coords; // What we're currently calculating
  calculation_metrics metrics;

  std::future<void> calculation_future;
  std::future<void> animation_future;
  std::atomic<bool> stop_calculation;
  std::atomic<bool> stop_animation;
  std::atomic<bool> calculation_completed;
  std::mutex mutex;
  std::condition_variable animation_condition;

  // True if we are calculating in the background.
  // This means that we'll display a zoom of the background image,
  // not the currently calculated image.
  bool animating;

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
};

// Perform a pixel-by-pixel remapping, interpolate values and don't increase the
// error values
void map_values(const view_pixmap &src, view_pixmap &dest, double dx, double dy,
                double r);

// Perform a pixel-by-pixel remapping, interpolate values and don't increase the
// error values
void interpolate_values(const view_pixmap &src, view_pixmap &dest, double dx,
                        double dy, double r);

} // namespace fractals
