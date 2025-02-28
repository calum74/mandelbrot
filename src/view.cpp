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
  if (!animating) {
    // We are rendering directly to the output
    // Copy the pixels to the output, update the metrics and notify the
    // listener.

    values = current_calculation_values;
    metrics.points_calculated = points_calculated;
    metrics.min_depth = min_depth;
    metrics.max_depth = max_depth;
    metrics.fully_evaluated = stride == 1;

    current_listener->update(metrics);
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

void fractals::view::animate_to(int x, int y,
                                std::chrono::duration<double> duration) {
  stop_calculating();
  stop_animating();

  previous_calculation_values = current_calculation_values;
}

void fractals::view::animate_to_center(std::chrono::duration<double> duration) {

}

void fractals::view::zoom(int x, int y, double r) {
  std::cout << "Zoom\n";
  freeze_current_view();

  current_coords.zoom(current_calculation_values.width(),
                      current_calculation_values.height(), x, y, r);
  previous_calculation_values = current_calculation_values;
  interpolate_values(previous_calculation_values, current_calculation_values, x,
                     y, r);
}

void fractals::view::scroll(int dx, int dy) {
  std::cout << "Scroll\n";
  freeze_current_view();

  current_coords.scroll(current_calculation_values.width(),
                        current_calculation_values.height(), dx, dy);

  // Scroll the pixels
  // In theory, we can move these within the same buffer
  // But this is lazier and we'll copy the buffer instead :-/

  previous_calculation_values = current_calculation_values;
  map_pixmap(
      previous_calculation_values, current_calculation_values, dx, dy, 1.0,
      [](auto c) { return c; }, missing_value);
}

void fractals::view::freeze_current_view() {
  // Stop everything and copy whatever we have into the current calculation
  stop_calculating();
  stop_animating();
  current_calculation_values = values;

  // TODO: Also remap the current coords.
}

namespace fractals {
template <typename ErrorFn>
void interpolate_values(const view_pixmap &src, view_pixmap &dest, double dx,
                        double dy, double r, ErrorFn fn) {

  for (int j = 0; j < dest.height(); ++j)
    for (int i = 0; i < dest.width(); ++i) {
      double rx = r * i + dx, ry = r * j + dy;
      int i2 = rx;
      int j2 = ry;
      auto &to_pixel = dest(i, j);
      if (i2 >= 0 && i2 < dest.width() && j2 >= 0 && j2 < dest.height()) {
        rx -= i2;
        ry -= j2;
        auto &from_pixel_00 = src(i2, j2);
        auto &from_pixel_10 = src(i2 + 1 < dest.width() ? i2 + 1 : i2, j2);
        auto &from_pixel_01 = src(i2, j2 + 1 < dest.height() ? j2 + 1 : j2);
        auto &from_pixel_11 = src(i2 + 1 < dest.width() ? i2 + 1 : i2,
                                  j2 + 1 < dest.height() ? j2 + 1 : j2);

        auto to_value = from_pixel_00.value * (1 - rx) * (1 - ry) +
                        from_pixel_10.value * rx * (1 - ry) +
                        from_pixel_01.value * (1 - rx) * ry +
                        from_pixel_11.value * rx * ry;
        auto to_error =
            std::max(std::max(from_pixel_00.error, from_pixel_10.error),
                     std::max(from_pixel_01.error, from_pixel_11.error));
        to_pixel = {to_value, fn(to_error)};
      } else {
        to_pixel = missing_value;
      }
    }
}
} // namespace fractals

void fractals::map_values(const view_pixmap &src, view_pixmap &dest, double dx,
                          double dy, double r) {
  bool zoom_eq = r == 1.0;
  bool zoom_out = r > 1.0;

  interpolate_values(src, dest, dx, dy, r, [&](int e) { return e; });
}

void fractals::interpolate_values(const view_pixmap &src, view_pixmap &dest,
                                  double dx, double dy, double r) {
  bool zoom_eq = r == 1.0;
  bool zoom_out = r > 1.0;

  interpolate_values(src, dest, dx, dy, r, [&](int e) {
    return zoom_eq ? e : zoom_out ? 20 : e > 20 ? e : e + 1;
  });
}
