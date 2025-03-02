#include "view.hpp"
#include "view_listener.hpp"

namespace fractals {
class view_animation : public view_listener {
public:
  view_animation(view_listener &listener);

  // Starts navigation at cursor, or aborts current animation
  void navigate_at_cursor(int x, int y);

  // Update the current position of the mouse
  void mouse_at(int x, int y);

  // Starts a single-step zoom at the cursor, or aborts the current animation
  void smooth_zoom_at_cursor(int x, int y);

  // Starts animating to the current position.
  void animate_to_current_position();

  // Sets the mode to use
  void set_quality_mode();
  void set_smooth_mode();
  void set_fast_mode();

private:
  fractals::view view;
  view_listener &listener;

  int mouse_x, mouse_y;

  enum class animation {
    none,
    navigate_at_point,
    navigate_to_point,
    navigate_randomly,
    single_zoom,
  } mode;
};
} // namespace fractals
