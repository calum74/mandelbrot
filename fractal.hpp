#pragma once

#include <memory>
#include <vector>

#include "ViewCoords.hpp"
#include "plane.hpp"

/* Basic infrastructure for rendering fractals. */

namespace fractals {

struct ViewCoords;
class Viewport;
class ColourMap;
class Renderer;

// A list of available fractals

std::unique_ptr<Renderer> make_renderer();

// Just a default colour map
std::unique_ptr<ColourMap> make_colourmap();

class PointwiseCalculation {
public:
  using ViewCoords = fractals::ViewCoords;
  virtual ~PointwiseCalculation() = default;

  virtual double calculate(int x, int y) const = 0;
};

class PointwiseFractal {
public:
  virtual const char *name() const = 0;
  virtual std::unique_ptr<PointwiseCalculation>
  create(const ViewCoords &c, int x, int y, std::atomic<bool> &stop) const = 0;
  virtual ViewCoords initial_coords() const = 0;
  virtual bool valid_for(const ViewCoords &c) const = 0;
};

namespace detail {

template <typename... Ts> class MultiPrecisionFactory;

template <typename T, typename... Ts>
class MultiPrecisionFactory<T, Ts...> : public PointwiseFractal {
public:
  MultiPrecisionFactory(const char *name) : n{name}, tail{name} {}

  ViewCoords initial_coords() const override { return T::initial_coords(); }
  const char *name() const override { return n; }

  std::unique_ptr<PointwiseCalculation>
  create(const ViewCoords &c, int x, int y,
         std::atomic<bool> &stop) const override {
    return T::valid_for(c) ? std::make_unique<T>(c, x, y, stop)
                           : tail.create(c, x, y, stop);
  }

  bool valid_for(const ViewCoords &c) const override {
    return T::valid_for(c) || tail.valid_for(c);
  }

private:
  MultiPrecisionFactory<Ts...> tail;

  const char *n;
};

template <typename T> class MultiPrecisionFactory<T> : public PointwiseFractal {
public:
  MultiPrecisionFactory(const char *name) : n{name} {}

  std::unique_ptr<PointwiseCalculation>
  create(const ViewCoords &c, int x, int y,
         std::atomic<bool> &stop) const override {
    return std::make_unique<T>(c, x, y, stop);
  }

  ViewCoords initial_coords() const override { return T::initial_coords(); }

  const char *name() const override { return n; }

  bool valid_for(const ViewCoords &c) const override { return T::valid_for(c); }

private:
  const char *n;
};
} // namespace detail

template <typename... Ts>
detail::MultiPrecisionFactory<Ts...> make_fractal(const char *name) {
  return {name};
}

class Registry {
public:
  virtual ~Registry() = default;
  virtual void add(const PointwiseFractal &) = 0;
  virtual std::vector<
      std::pair<std::string, const fractals::PointwiseFractal &>>
  listFractals() const = 0;
};

// Perhaps this isn't useful
std::unique_ptr<Registry> make_registry();

}; // namespace fractals

