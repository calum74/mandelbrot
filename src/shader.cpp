#include "shader.hpp"
#include <cmath>
#include <random>

double fractals::calculate_brightness(double dx, double dy,
                                      double colour_gradient,
                                      double ambient_brightness,
                                      double source_brightness,
                                      const unit_vector &light_source) {
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

std::vector<fractals::RGB> fractals::generate_colours(int size, int seed) {
  std::mt19937 e(seed);
  std::vector<RGB> newColours(size);
  for (auto &c : newColours) {
    c = e() & 0xffffff;
  }
  return newColours;
}

fractals::unit_vector fractals::spherical_to_cartesian(double direction,
                                                       double elevation) {
  return {std::cos(elevation) * std::cos(direction),
          std::cos(elevation) * std::sin(direction), std::sin(elevation)};
}

fractals::RGB fractals::get_colour_from_index(const std::vector<RGB> &colours,
                                              double index, double brightness) {
  int i = index;
  auto f2 = index - i;
  auto f1 = 1.0 - f2;

  i %= colours.size();
  int j = (i + 1) % colours.size();
  auto c1 = colours[i];
  auto c2 = colours[j];

  auto r = brightness * (red(c1) * f1 + red(c2) * f2);
  auto g = brightness * (green(c1) * f1 + green(c2) * f2);
  auto b = brightness * (blue(c1) * f1 + blue(c2) * f2);

  if (r > 255)
    r = 255;
  if (g > 255)
    g = 255;
  if (b > 255)
    b = 255;

  if (r < 0)
    r = 0;
  if (g < 0)
    g = 0;
  if (b < 0)
    b = 0;

  return make_rgb(r, g, b);
}
