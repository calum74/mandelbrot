#include "fractal.hpp"
#include <string>
#include <vector>

namespace {
using namespace fractals;

class RegistryImpl : public Registry {
  void add(const PointwiseFractal &f) override {
    fractals.push_back(
        std::pair<std::string, const PointwiseFractal &>{f.name(), f});
  }

  std::vector<std::pair<std::string, const fractals::PointwiseFractal &>>
      fractals;

  std::vector<std::pair<std::string, const fractals::PointwiseFractal &>>
  listFractals() const override {
    return fractals;
  }
};
} // namespace

std::unique_ptr<fractals::Registry> fractals::make_registry() {
  return std::make_unique<RegistryImpl>();
}
