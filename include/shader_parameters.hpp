#pragma once

namespace fractals {

/*
    Options that control the rendering of the fractal.
*/
struct shader_parameters {
  int colour_scheme =168; // Any integer used to seed the random number generator
  double colour_gradient = 30; // Iterations per colour
  double colour_offset = 0;
  bool auto_gradient = true;

  bool shading = true;
  double ambient_brightness = 0.4;       // Value between 0 and 1
  double source_brightness = 0.5;        // Value between 0 and 1
  double source_direction_radians = 3.14/4; // Value between 0 and 2π
  double source_elevation_radians = 3.14/4; // Value between 0 and π
};
} // namespace fractals
