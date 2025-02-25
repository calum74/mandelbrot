#pragma once

namespace fractals {

/*
    Options that control the rendering of the fractal.
*/
struct shader_parameters {
  int colour_scheme; // Any integer used to seed the random number generator
  double colour_gradient; // Iterations per colour
  double colour_offset;
  bool auto_gradient;

  bool shading;
  double ambient_brightness;       // Value between 0 and 1
  double source_brightness;        // Value between 0 and 1
  double source_direction_radians; // Value between 0 and 2π
  double source_elevation_radians; // Value between 0 and π
};
} // namespace fractals
