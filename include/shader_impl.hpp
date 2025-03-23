#include "shader.hpp"
#include "shader_parameters.hpp"
#include "gradient_stack.hpp"

namespace fractals {

class shader_impl : public shader {
public:
  shader_impl();
  RGB operator()(double d, double dx, double dy) const override;
  RGB operator()(double d) const override;
  void randomize() override;
  void resetGradient() override;
  void setRange(double min, double max) override;
  void maybeUpdateRange(double min, double max) override;
  void load(const view_parameters &) override;
  void save(view_parameters &) const override;

  void setParameters(const shader_parameters &) override;
  void getParameters(shader_parameters &params) override;

private:
  void create_colours();
  void update_light_source();

  shader_parameters params;
  const unsigned int numColours = 100;
  std::vector<RGB> colours;
  gradient_stack gradients;
  unit_vector light_source;
};
} // namespace fractals
