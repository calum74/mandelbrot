#include "shader.hpp"
#include <cmath>
#include <random>

double fractals::calculate_brightness(double dx, double dy,
                                      double colour_gradient,
                                      double ambient_brightness,
                                      double source_brightness, const unit_vector & light_source) {
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
      (surface_normal_x * light_source.x + surface_normal_y * light_source.y +
       surface_normal_z * light_source.z) /
      surface_normal_length;

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
