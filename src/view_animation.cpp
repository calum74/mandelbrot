#include "view_animation.hpp"

using namespace std::literals::chrono_literals;

fractals::view_animation::view_animation(fractals::view_listener &listener)
    : listener(listener) {
  view.set_listener(this);
  wait_for_completion = false;
  zoom_step_duration = 250ms;
  navigate_step_duration = 750ms;
  animate_step_duration = 500ms;

  mode = animation::none;
  animation_loop_thread = std::async([&] { animation_loop(); });
}

fractals::view_animation::~view_animation() {
  mode = animation::shutdown;
  animation_variable.notify_one();
  animation_loop_thread.wait();
}

void fractals::view_animation::animate_to(
    int x, int y, std::chrono::duration<double> duration,
    bool wait_for_completion) {
  view.animate_to(x, y, duration, wait_for_completion);
}

void fractals::view_animation::set_coords(const view_coords &coords,
                                          bool recalculate) {
  view.set_coords(coords, recalculate);
}

void fractals::view_animation::set_fractal(const fractal &f, bool init,
                                           bool recalculate) {
  view.set_fractal(f, init, recalculate);
}

void fractals::view_animation::set_threading(int threading) {
  view.set_threading(threading);
}

void fractals::view_animation::animate_to_center(
    std::chrono::duration<double> duration, bool wait_for_completion) {
  view.animate_to_center(duration, wait_for_completion, 0.5);
}

void fractals::view_animation::start_calculating() { view.start_calculating(); }

void fractals::view_animation::update_iterations(
    fractals::calculation_metrics const &metrics) {
  view.update_iterations(metrics);
}

void fractals::view_animation::set_max_iterations(int max) {
  view.set_max_iterations(max);
}

void fractals::view_animation::decrease_iterations() {
  view.decrease_iterations();
}

void fractals::view_animation::increase_iterations() {
  view.increase_iterations();
}

void fractals::view_animation::stop_current_animation_and_set_as_current() {
  mode = animation::none;
  view.stop_current_animation_and_set_as_current();
}

void fractals::view_animation::zoom(int x, int y, double ratio) {
  view.zoom(x, y, ratio);
}

void fractals::view_animation::scroll(int dx, int dy) { view.scroll(dx, dy); }

void fractals::view_animation::set_size(int w, int h) {
  // TODO: Could be delayed until the next calculation
  view.set_size(w, h);
}

const fractals::view_coords &fractals::view_animation::get_coords() const {
  return view.get_coords();
}

const fractals::calculation_metrics &
fractals::view_animation::get_metrics() const {
  return view.get_metrics();
}

const fractals::view_pixmap &fractals::view_animation::values() const {
  return view.values;
}

bool fractals::view_animation::is_animating() const {
  return view.is_animating();
}

fractals::view_coords fractals::view_animation::initial_coords() const {
  return view.initial_coords();
}

std::string fractals::view_animation::get_fractal_name() const {
  return view.get_fractal_name();
}

std::string fractals::view_animation::get_fractal_family() const {
  return view.get_fractal_family();
}

int fractals::view_animation::width() const { return view.width(); }

int fractals::view_animation::height() const { return view.height(); }

void fractals::view_animation::calculation_started(radius r,
                                                   int max_iterations) {
  listener.calculation_started(r, max_iterations);
}

void fractals::view_animation::values_changed() { listener.values_changed(); }

void fractals::view_animation::calculation_finished(
    const calculation_metrics &metrics) {
  listener.calculation_finished(metrics);

  if (mode == animation::start_navigate_to_point) {
    mode = animation::navigate_to_point;
    animation_variable.notify_one();
  }
}

void fractals::view_animation::animation_finished(
    const calculation_metrics &metrics) {
  listener.animation_finished(metrics);
  animation_variable.notify_one();
}

void fractals::view_animation::navigate_at_cursor(int x, int y) {
  std::unique_lock<std::mutex> lock(mutex);

  mouse_at(x, y);

  if (mode == animation::none) {
    std::cout << ":: Starting navigate\n";
    mode = animation::navigate_at_cursor;
    view.animate_to(x, y, navigate_step_duration, wait_for_completion);
  } else {
    std::cout << ":: Stopping navigate\n";
    mode = animation::none;
    view.stop_current_animation_and_set_as_current();
  }
}

void fractals::view_animation::mouse_at(int x, int y) {
  mouse_x = x;
  mouse_y = y;
}

void fractals::view_animation::animation_loop() {
  std::unique_lock<std::mutex> lock(mutex);
  while (mode != animation::shutdown) {
    animation_variable.wait(lock);
    std::cout << "Animation loop notification received\n";
    switch (mode) {
    case animation::navigate_at_cursor:
      view.animate_to(mouse_x, mouse_y, navigate_step_duration,
                      wait_for_completion);
      break;
    case animation::start_navigate_to_point:
      mode = animation::navigate_to_point;
      // Fall through
    case animation::navigate_to_point: {
      auto current_r = view.get_coords().ln_r();
      auto new_r = current_r - std::log(2);
      if (new_r > zoom_limit)
        view.animate_to_center(animate_step_duration, wait_for_completion, 0.5);
      else {
        auto r = std::exp(zoom_limit - current_r);
        // TODO: Adjust animation duration on the last step
        view.animate_to_center(animate_step_duration, wait_for_completion, r);
        mode = animation::none;
      }
    } break;
    case animation::none:
    case animation::navigate_randomly:
    case animation::single_zoom:
    case animation::shutdown:
      // TODO
      break;
    }
  }
}

void fractals::view_animation::animate_to_current_position() {
  std::unique_lock<std::mutex> lock(mutex);
  mode = animation::start_navigate_to_point;

  auto c = view.get_coords();
  zoom_limit = c.ln_r();
  c.r = 2.0;
  c.max_iterations = 500;
  view.set_coords(c, true);
}
