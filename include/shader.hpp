#pragma once
#include "rgb.hpp"
#include <vector>

namespace fractals {
double calculate_brightness(double dx, double dy, double colour_gradient,
                            double ambient_brightness, double source_brightness,
                            double source_x, double source_y, double source_z,
                            double source_length);

std::vector<RGB> generate_colours(int size, int seed);

} // namespace fractals