#include "shader.hpp"
#include <cmath>
#include <random>

double fractals::calculate_brightness(double dx, double dy,
                                      double colour_gradient,
                                      double ambient_brightness,
                                      double source_brightness, double source_x,
                                      double source_y, double source_z,
                                      double source_length) {
  dx /= colour_gradient;
  dy /= colour_gradient;

  dx *= 3000;
  dy *= 3000;

  double surface_normal_x = -dx;
  double surface_normal_y = -dy;
  double surface_normal_z = 1;
  double surface_normal_length = std::sqrt(surface_normal_x * surface_normal_x +
                                           surface_normal_y * surface_normal_y +
                                           surface_normal_z * surface_normal_z);

  double dot_product =
      (surface_normal_x * source_x + surface_normal_y * source_y +
       surface_normal_z * source_z) /
      (surface_normal_length * source_length);

  if (dot_product < 0)
    dot_product = 0;

  return ambient_brightness + source_brightness * dot_product;
}

std::vector<fractals::RGB> fractals::generate_colours(int size, int seed)
{
  std::mt19937 e(seed);
  std::vector<RGB> newColours(size);
  for (auto &c : newColours) {
    c = e() & 0xffffff;
  }
  return newColours;
}
