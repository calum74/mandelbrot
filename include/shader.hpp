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

struct view_parameters;
struct shader_parameters;

class shader {
public:
  virtual ~shader() = default;
  virtual void randomize() = 0;
  virtual void resetGradient() = 0;
  virtual void setRange(double min, double max) = 0;
  virtual void maybeUpdateRange(double min, double max) = 0;
  virtual RGB operator()(double depth, double dx, double dy) const = 0;
  virtual RGB operator()(double depth) const = 0;

  virtual void load(const view_parameters &) = 0;
  virtual void save(view_parameters &) const = 0;
  virtual void setParameters(const shader_parameters &) = 0;
  virtual void getParameters(shader_parameters &params) = 0;
};

std::unique_ptr<shader> make_shader();

} // namespace fractals