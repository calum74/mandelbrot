#pragma once
#include "rgb.hpp"
#include <vector>

namespace fractals {

struct unit_vector {
  double x, y, z;
};

unit_vector spherical_to_cartesian(double direction_radians, double incline_radians);

double calculate_brightness(double dx, double dy, double colour_gradient,
                            double ambient_brightness, double source_brightness,
                            const unit_vector & light_source);

std::vector<RGB> generate_colours(int size, int seed);

RGB get_colour_from_index(const std::vector<RGB> &colours, double index, double brightness);

} // namespace fractals