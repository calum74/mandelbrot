#include "shader_impl.hpp"
#include "view_parameters.hpp"

fractals::shader_impl::shader_impl() {
  resetGradient();
  create_colours();
  update_light_source();
}

void fractals::shader_impl::resetGradient() {
  params.colour_gradient = 30;
  params.colour_offset = 0;
  gradients.clear();
}

fractals::RGB fractals::shader_impl::operator()(double d, double dx,
                                                  double dy) const {
  if (d == 0)
    return make_rgb(0, 0, 0);

  if (!params.shading)
    return (*this)(d);

  auto colour_index =
      gradients.map_iteration(d, params.colour_gradient, params.colour_offset);

  double brightness = calculate_brightness(
      dx, dy, colour_index.gradient, params.ambient_brightness,
      params.source_brightness, light_source);

  return get_colour_from_index(colours, colour_index.value, brightness);
}

fractals::RGB fractals::shader_impl::operator()(double d) const {
  if (d == 0)
    return make_rgb(0, 0, 0);

  auto colour_index =
      gradients.map_iteration(d, params.colour_gradient, params.colour_offset);

  return get_colour_from_index(colours, colour_index.value, 1.0);
}

void fractals::shader_impl::setRange(double min, double max) {
  params.colour_gradient = (max - min) / 5.0;
  gradients.clear();
}

void fractals::shader_impl::maybeUpdateRange(double min, double max) {

  if (!params.auto_gradient)
    return;

  gradients.push(max, (max - min) / 5, params.colour_gradient,
                 params.colour_offset);
}

void fractals::shader_impl::randomize() {
  params.colour_scheme++;
  create_colours();
}

void fractals::shader_impl::create_colours() {
  colours = generate_colours(numColours, params.colour_scheme);
}

void fractals::shader_impl::load(const view_parameters &vp) {
  params.colour_scheme = vp.shader.colour_scheme;
  params.colour_gradient = vp.shader.colour_gradient;
  if (params.colour_gradient < 1.0)
    params.colour_gradient = 1.0 / params.colour_gradient;
  gradients.clear();
  create_colours();
}

void fractals::shader_impl::save(view_parameters &vp) const {
  vp.shader = params;
}

void fractals::shader_impl::getParameters(shader_parameters &vp) {
  vp = params;
}

std::unique_ptr<fractals::shader> fractals::make_shader() {
  return std::make_unique<shader_impl>();
}

void fractals::shader_impl::setParameters(const shader_parameters &vp) {
  if (params.colour_scheme != vp.colour_scheme) {
    params.colour_scheme = vp.colour_scheme;
    create_colours();
  }

  if (params.colour_gradient != vp.colour_gradient ||
      params.colour_offset != vp.colour_offset) {
    gradients.clear();
  }

  if (params.auto_gradient && !vp.auto_gradient) {
    gradients.clear();
  }

  params = vp;
  update_light_source();
}

void fractals::shader_impl::update_light_source() {
  // Recalculate the source based on spherical coordinates
  light_source = spherical_to_cartesian(params.source_direction_radians,
                                        params.source_elevation_radians);
}
