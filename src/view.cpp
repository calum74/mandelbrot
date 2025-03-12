#include "view.hpp"
#include "calculation_pixmap.hpp"
#include "fractal_calculation.hpp"
#include "percentile.hpp"
#include "view_listener.hpp"

using namespace std::literals::chrono_literals;

fractals::view::view() : listener{}, animating(false) {
  calculation_threads = 4;
}

fractals::view::~view() {
  stop_calculating();
  stop_animating();
}

void fractals::view::set_listener(view_listener *l) { listener = l; }

void fractals::view::stop_calculating() {
  if (calculation_future.valid()) {
    stop_calculation = true;
    calculation_future.wait();
    stop_calculation = false;
  }
}

struct fractals::view::my_calculation_pixmap
    : public fractals::calculation_pixmap {
  fractals::view &view;

  my_calculation_pixmap(fractals::view &v)
      : fractals::calculation_pixmap(v.current_calculation_values, 16,
                                     *v.calculation),
        view(v) {}

  virtual void layer_complete(int stride) {
    view.complete_layer(min_depth, max_depth, points_calculated, stride);
  }
};

void fractals::view::start_calculating() {
  if (!valid())
    return;
  stop_calculating();
  calculation_future = std::async([&] {
    auto start_time = std::chrono::system_clock::now();
    {
      std::unique_lock<std::mutex> m(mutex);
      metrics.log_radius = calculation_coords.ln_r();
      std::cout << "  Starting calculating with ln_r = " << metrics.log_radius
                << std::endl;
      calculation_completed = false;
      listener->calculation_started(metrics.log_radius,
                                    calculation_coords.max_iterations);
    }

    calculation = fractal->create(calculation_coords, width(), height(),
                                  stop_calculation);

    my_calculation_pixmap pm(*this);
    pm.calculate(calculation_threads, stop_calculation);

    metrics.fully_evaluated = calculation_completed;

    std::chrono::duration<double> duration =
        std::chrono::system_clock::now() - start_time;

    metrics.seconds_per_point = duration.count() / metrics.points_calculated;
    metrics.render_time_seconds = duration.count();
    metrics.discovered_depth = metrics.max_depth;
    measure_depths(current_calculation_values, metrics);
    listener->calculation_finished(metrics);
    std::cout << "  Completed = " << calculation_completed << "\n";
    std::cout << "==== End of calculation\n\n";
  });
}

void fractals::view::stop_animating() {
  if (animation_future.valid()) {
    stop_animation = true;
    animation_condition.notify_one();
    animation_future.wait();
    stop_animation = false;
  }
  animating = false;
}

void fractals::view::animation_thread() {
  while (!stop_animation) {
    std::unique_lock<std::mutex> lock(mutex);
    auto e = animation_condition.wait_for(lock, 10ms);

    if (e == std::cv_status::timeout && !stop_animation) {
      std::chrono::duration<double> duration =
          std::chrono::system_clock::now() - animation_start;

      double time_ratio = duration.count() / animation_duration.count();

      if (time_ratio >= 1) {
        time_ratio = 1;
        stop_animation = true;
        animating = false;
      }
      rendered_zoom_ratio = std::pow(current_step_ratio, time_ratio);

      // !! Use the mutex rather than atomics
      bool display_previous_image =
          time_ratio < 1 ||
          (wait_for_calculation_to_complete && !calculation_completed);

      if (display_previous_image) {

        map_values(previous_calculation_values, values,
                   zoom_x * (1.0 - rendered_zoom_ratio),
                   zoom_y * (1.0 - rendered_zoom_ratio), rendered_zoom_ratio);

      } else {
        std::cout << "Animation finished\n";

        // Animation finished
        if (calculation_completed) {
          std::cout << "Calculation completed\n";

          values = current_calculation_values;

        } else {
          std::cout << "Waiting for calculation\n";
        }
        // Except in quality mode
        listener->animation_finished(metrics); // Maybe start another animation
        // else: Wait for the calculation to update the image
      }
      listener->values_changed();
    }
  }
  std::cout << "Finished animating\n";
}

void fractals::view::complete_layer(double min_depth, double max_depth,
                                    std::uint64_t points_calculated,
                                    int stride) {
  std::cout << "  Layer completed\n";
  std::unique_lock<std::mutex> lock(mutex);

  if (stride == 1) {
    calculation_completed = true;
    std::cout << "  calculation_completed = " << calculation_completed
              << std::endl;
  }

  if (!animating) {
    // We are rendering directly to the output
    // Copy the pixels to the output, update the metrics and notify the
    // listener.

    values = current_calculation_values;
    metrics.points_calculated = points_calculated;
    std::cout << "  " << points_calculated << " points calculated\n";
    metrics.min_depth = min_depth;
    metrics.max_depth = max_depth;
    metrics.fully_evaluated = stride == 1;
    listener->values_changed();
  }
}

constexpr fractals::error_value<double> missing_value{
    std::numeric_limits<double>::quiet_NaN(), 127};

void fractals::view::set_size(int w, int h) {
  stop_animating();
  stop_calculating();

  // TODO: Recalculate aspect ratio and coords

  // TODO: Remap the viewport smoothly
  values = view_pixmap(w, h, missing_value);
  previous_calculation_values = view_pixmap(w, h, missing_value);
  current_calculation_values = view_pixmap(w, h, missing_value);
}

void fractals::view::set_threading(int n) {
  calculation_threads = n > 0 ? n : std::thread::hardware_concurrency();
}

void fractals::view::set_fractal(const fractals::fractal &f, bool init_coords,
                                 bool recalculate) {
  if (!init_coords && f.name() == get_fractal_name())
    return;

  fractal = f.create();
  if (init_coords)
    set_coords(fractal->initial_coords(), recalculate);
  else if (recalculate)
    start_calculating();
}

void fractals::view::set_coords(const view_coords &vc, bool recalculate) {
  stop_animating();
  stop_calculating();
  calculation_coords = vc;
  invalidate_values(current_calculation_values);
  if (recalculate)
    start_calculating();
}

bool fractals::view::valid() const {
  return fractal && listener && calculation_threads > 0 &&
         current_calculation_values.width() > 0 &&
         current_calculation_values.height() > 0;
}

void fractals::view::animate_to(int x, int y,
                                std::chrono::duration<double> duration,
                                bool wait, bool lock_center, double ratio) {
  stop_calculating();
  stop_animating();

  previous_calculation_values = current_calculation_values;

  animating = true;
  rendered_zoom_ratio = 1;
  animation_start = std::chrono::system_clock::now();
  animation_duration = duration;
  calculation_completed = false;
  std::cout << "  calculation_completed = " << calculation_completed
            << std::endl;
  wait_for_calculation_to_complete = wait;
  zoom_x = x;
  zoom_y = y;

  // Where to calculate
  current_step_ratio = ratio;
  double r = ratio;
  calculation_coords = lock_center ? calculation_coords.zoom(r) : calculation_coords.zoom(r, width(), height(), x, y);

  // Seed the current calculation values so that if we abort, we already have
  // some data there already
  map_values(values, current_calculation_values, x * (1 - r), y * (1 - r), r);

  // Start the animation thread
  stop_animation = false;
  animation_future = std::async([&]() { animation_thread(); });

  // Start the calculation thread
  start_calculating();
}

void fractals::view::animate_to(int x, int y,
                                std::chrono::duration<double> duration,
                                bool wait) {
  animate_to(x, y, duration, wait, false, 0.5);
}

void fractals::view::animate_to_center(std::chrono::duration<double> duration,
                                       bool wait, double ratio) {
  animate_to(width() / 2, height() / 2, duration, wait, true, ratio);
}

void fractals::view::zoom(int x, int y, double r) {
  std::cout << "Zoom\n";

  if (r == 1.0)
    return;
  if (r > 2)
    r = 2;
  if (r < 0.5)
    r = 0.5;

  freeze_current_view();

  calculation_coords = calculation_coords.zoom(r, width(), height(), x, y);

  map_values(values, current_calculation_values, x * (1 - r), y * (1 - r), r);
  values = current_calculation_values;
  listener->values_changed();
  start_calculating();
}

void fractals::view::scroll(int dx, int dy) {

  if (dx == 0 && dy == 0)
    return;

  freeze_current_view();

  calculation_coords = calculation_coords.scroll(width(), height(), dx, dy);

  // Scroll the pixels
  // In theory, we can move these within the same buffer
  // But this is lazier and we'll copy the buffer instead :-/

  previous_calculation_values = current_calculation_values;
  map_pixmap(
      previous_calculation_values, current_calculation_values, dx, dy, 1.0,
      [](auto c) { return c; }, missing_value);
  values = current_calculation_values;
  listener->values_changed();
  start_calculating();
}

void fractals::view::freeze_current_view() {
  // Stop everything and copy whatever we have into the current calculation
  stop_calculating();
  stop_animating();

  if (animating) {
    std::cout << "Animating??\n";
    current_calculation_values = values;
    // TODO: Also remap the current coords.
  }
}

bool fractals::view::get_auto_zoom(int &, int &) {
  // TODO
  return false;
}

void fractals::view::update_iterations(const calculation_metrics &) {
  // TODO
}

void fractals::view::decrease_iterations() {
  // TODO
}

void fractals::view::increase_iterations() {
  // TODO
}

void fractals::view::stop_current_animation_and_set_as_current() {
  std::cout << "Aborting animation\n";
  stop_calculating();
  stop_animating();

  // TODO: Update coords and copy.

  calculation_coords = calculation_coords.zoom(
      (1.0/current_step_ratio) * rendered_zoom_ratio, width(), height(), zoom_x, zoom_y);
  current_calculation_values = values;
  start_calculating();
}

const fractals::view_coords &fractals::view::get_coords() const {
  return calculation_coords;
}

const fractals::calculation_metrics &fractals::view::get_metrics() const {
  return metrics;
}

fractals::view_coords fractals::view::initial_coords() const {
  return fractal->initial_coords();
}

std::string fractals::view::get_fractal_name() const { return fractal->name(); }

std::string fractals::view::get_fractal_family() const {
  return fractal->family();
}

void fractals::view::set_max_iterations(int max) {
  calculation_coords.max_iterations = max;
}

void fractals::measure_depths(const view_pixmap &values,
                              calculation_metrics &metrics) {
  // Find the depths.
  std::vector<int> coloured_pixels;
  coloured_pixels.reserve(values.size());
  for (int i = 0; i < values.size(); ++i) {
    if (values[i].error == 0 && values[i].value > 0)
      coloured_pixels.push_back(i);
  }

  metrics.non_black_points = coloured_pixels.size();

#if 0
  std::cout << "Finished calculation\n";
  std::cout << "  We calculated " << calculated_pixels.points_calculated
            << " points\n";

  std::cout << "  We have " << coloured_pixels.size() << " coloured pixels\n";
#endif

  if (coloured_pixels.size() > 100) {
    // ?? Why are we doing this? Can't we just take the top pixel and be done
    // with it?
    auto cmp = [&](int a, int b) { return values[a].value < values[b].value; };
    auto p999 = *top_percentile(coloured_pixels.begin(), coloured_pixels.end(),
                                0.999, cmp);
    auto p9999 = *top_percentile(coloured_pixels.begin(), coloured_pixels.end(),
                                 0.9999, cmp);
    metrics.p999 = values[p999].value;
    metrics.p9999 = values[p9999].value;
    metrics.p9999_x = p9999 % values.width();
    metrics.p9999_y = p9999 / values.width();
  }
}

bool fractals::view::is_animating() const { return animating; }