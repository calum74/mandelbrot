#include "view.hpp"
#include "view_listener.hpp"
#include "fractal_calculation.hpp"

namespace fractals {

class view_animation : public view_listener {
public:
  view_animation(view_listener &listener);
  ~view_animation();

  // Starts navigation at cursor, or aborts current animation
  void navigate_at_cursor(int x, int y);

  // Update the current position of the mouse
  void mouse_at(int x, int y);

  // Starts a single-step zoom at the cursor, or aborts the current animation
  void smooth_zoom_at_cursor(int x, int y);

  void scroll(int dx, int dy);

  void zoom(int x, int y, double ratio);

  // Starts animating to the current position.
  void animate_to_current_position();

  void navigate_randomly();

  // Sets the mode to use
  void set_quality_mode();
  void set_smooth_mode();
  void set_fast_mode();

  void set_threading(int);
  void set_fractal(const fractal &, bool, bool);
  int width() const;
  int height() const;

  void animate_to(int x, int y, std::chrono::duration<double> duration,
                  bool wait_for_completion);
  const view_coords &get_coords() const;
  void set_coords(const view_coords &, bool recalculate);

  void set_max_iterations(int);
  void increase_iterations();
  void decrease_iterations();
  const view_pixmap &values() const;

  std::string get_fractal_name() const;
  std::string get_fractal_family() const;

  void set_size(int w, int h);

  view_coords initial_coords() const;
  const calculation_metrics &get_metrics() const;

  void get_orbit(int x, int y, displayed_orbit &orbit) const;

  // Deleteme!
  void animate_to_center(std::chrono::duration<double> duration,
                         bool wait_for_completion);
  bool is_animating() const;
  void stop_current_animation_and_set_as_current();
  void update_iterations(const calculation_metrics &);
  void start_calculating();

  // If the current values (interpolated or not) are a result of a completed
  // calculation
  bool fully_calculated() const;

  // Can just assign to these directly as needed:
  std::chrono::duration<double> zoom_step_duration, navigate_step_duration,
      animate_step_duration;
  bool wait_for_completion;

private:
  fractals::view view;
  view_listener &listener;

  std::future<void> animation_loop_thread;
  std::mutex mutex;
  std::condition_variable animation_variable;
  radius zoom_limit;

  std::chrono::duration<double>
      quality_duration; // The best guess how long it takes to calculate a
                        // single frame

  std::pair<std::chrono::duration<double>, bool>
  get_step_duration(std::chrono::duration<double> d) const;

  // Protected by mutex
  int mouse_x, mouse_y;
  int random_x, random_y;

  // Protected by mutex
  enum class animation {
    none,
    navigate_at_cursor,
    start_navigate_to_point,
    navigate_to_point,
    navigate_randomly,
    single_zoom,
    shutdown
  } mode;

  void calculation_started(radius, int max_iterations) override;
  void values_changed() override;
  void calculation_finished(const calculation_metrics &) override;
  void animation_finished(const calculation_metrics &) override;

  void animation_loop();
};
} // namespace fractals
