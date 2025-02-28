#include "view.hpp"
#include "calculation_pixmap.hpp"
#include "fractal_calculation.hpp"
#include <chrono>

using namespace std::literals::chrono_literals;

fractals::view::view() : current_listener{} { calculation_threads = 4; }

fractals::view::~view() {
  stop_calculating();
  stop_animating();
}

void fractals::view::set_listener(listener *l) { current_listener = l; }

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
    my_calculation_pixmap pm(*this);
    pm.calculate(calculation_threads, stop_calculation);
  });
}

void fractals::view::stop_animating() {
  if (animation_future.valid()) {
    stop_animation = true;
    animation_condition.notify_one();
    animation_future.wait();
    stop_animation = false;
  }
}

void fractals::view::animation_thread() {
  while (!stop_animation) {
    std::unique_lock<std::mutex> lock;
    auto e = animation_condition.wait_for(lock, 10ms);

    if (e == std::cv_status::timeout && !stop_animation) {

      // If we are in the animation time,
      // Are we in the animation time
      // Look at the wallclock time...
      // Time to
      // The condition

      // If we have finished animating, then either
      // a) we have finished calculating
      // b) we have not finished calculating, and we need to wait for the
      // calculation to finish c) we haven't finished calculating, and we show
      // the intermediate results straight away. If we finished calculating,
      // then copy the pixels to the output and we are finished If we haven't
      // finished calculating, then
    }
  }
}

void fractals::view::complete_layer(double min_depth, double max_depth,
                                    std::uint64_t points_calculated,
                                    int stride) {
  if(!animating)
  {
    // We are rendering directly to the output
    // Copy the pixels to the output, update the metrics and notify the listener.

    for(std::size_t i=0; i<current_calculation_values.size(); ++i)
      values[i] = current_calculation_values[i].value;

      metrics.points_calculated = points_calculated;
    metrics.min_depth = min_depth;
    metrics.max_depth = max_depth;
    metrics.fully_evaluated = stride==1;

    current_listener->update(metrics);
  }
}

void fractals::view::set_size(int w, int h) {
  stop_animating();
  stop_calculating();

  // TODO: Recalculate aspect ratio and coords

  // TODO: Remap the viewport smoothly
  values = pixmap<double>(w, h, 0);
  error_value<double> missing_value{std::numeric_limits<double>::quiet_NaN(),
                                    127};
  previous_calculation_values =
      pixmap<error_value<double>>(w, h, missing_value);
  current_calculation_values = pixmap<error_value<double>>(w, h, missing_value);
}

void fractals::view::set_threading(int n) {
  calculation_threads = n > 0 ? n : std::thread::hardware_concurrency();
}

void fractals::view::set_fractal(const fractals::fractal &f) {
  fractal = f.create();
  set_coords(fractal->initial_coords());
}

void fractals::view::set_coords(const view_coords &vc) {
  stop_animating();
  stop_calculating();
  current_coords = fractal->initial_coords();
  calculation =
      fractal->create(current_coords, current_calculation_values.width(),
                      current_calculation_values.height(), stop_calculation);
}

bool fractals::view::valid() const {
  return fractal && current_listener && calculation_threads > 0 &&
         current_calculation_values.width() > 0 &&
         current_calculation_values.height() > 0;
}

void fractals::view::animate_to(int x, int y, std::chrono::duration<double> duration)
{
  stop_calculating();
  stop_animating();

  previous_calculation_values = current_calculation_values;
}

void fractals::view::animate_to_center(std::chrono::duration<double> duration)
{
  
}

void fractals::view::zoom(int x, int y, double r)
{
  std::cout << "Zoom\n";
}

void fractals::view::scroll(int dx, int dy)
{
  std::cout << "Scroll\n";
}