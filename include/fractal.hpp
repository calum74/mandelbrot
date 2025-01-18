#pragma once

#include <memory>

#include "plane.hpp"
#include "view_coords.hpp"

namespace fractals {

/*
  The calculation of a fractal, where we can compute each (x,y) point
  independently and perhaps in parallel.

  A virtual function call per point may seem excessive, but it's usually
  negligible overhead relative to computing each point.
*/
class PointwiseCalculation {
public:
  using view_coords = fractals::view_coords;
  virtual ~PointwiseCalculation() = default;

  // Compute the point at the pixel position (x,y).
  virtual double calculate(int x, int y) const = 0;

  // Gets the average number of iterations, for fractals where this makes sense.
  virtual double average_iterations() const;

  // Gets the average number of skipped iterations, for fractals where this
  // makes sense.
  virtual double average_skipped() const;

  virtual void initialize(const view_coords &c, int x, int y,
                          std::atomic<bool> &stop);
};

/*
  A fractal that can be calculated in a point-wise fashion.
*/
class PointwiseFractal {
public:
  virtual const char *name() const = 0;

  virtual const char *family() const = 0; //

  // Creates and initializes a new calculation.
  // Called once for each generated fractal image.
  virtual std::shared_ptr<PointwiseCalculation>
  create(const view_coords &c, int x, int y, std::atomic<bool> &stop) const = 0;

  // Retrieve the initial coordinates for this fractal.
  virtual view_coords initial_coords() const = 0;

  // Queries whether the given coordinates are valid for this fractal type.
  virtual bool valid_for(const view_coords &c) const = 0;
};

namespace detail {

template <typename... Ts> class MultiPrecisionFactory;

template <typename T, typename... Ts>
class MultiPrecisionFactory<T, Ts...> : public PointwiseFractal {
public:
  MultiPrecisionFactory(const char *name, const char *family)
      : n{name}, f{family}, tail{name, family} {}

  view_coords initial_coords() const override { return T::initial_coords(); }
  const char *name() const override { return n; }

  const char *family() const override { return f; }

  mutable std::shared_ptr<PointwiseCalculation> previous;

  std::shared_ptr<PointwiseCalculation>
  create(const view_coords &c, int x, int y,
         std::atomic<bool> &stop) const override {

    if (T::valid_for(c)) {
      if (!previous)
        previous = std::make_shared<T>(c, x, y, stop);
      previous->initialize(c, x, y, stop);
      return previous;
    }
    return tail.create(c, x, y, stop);
  }

  bool valid_for(const view_coords &c) const override {
    return T::valid_for(c) || tail.valid_for(c);
  }

private:
  MultiPrecisionFactory<Ts...> tail;

  const char *n;
  const char *f;
};

template <typename T> class MultiPrecisionFactory<T> : public PointwiseFractal {
public:
  MultiPrecisionFactory(const char *name, const char *family)
      : n{name}, f{family} {}

  mutable std::shared_ptr<PointwiseCalculation> previous;

  std::shared_ptr<PointwiseCalculation>
  create(const view_coords &c, int x, int y,
         std::atomic<bool> &stop) const override {
    if (!previous)
      previous = std::make_shared<T>(c, x, y, stop);
    previous->initialize(c, x, y, stop);
    return previous;
  }

  view_coords initial_coords() const override { return T::initial_coords(); }

  const char *name() const override { return n; }

  const char *family() const override { return f; }

  bool valid_for(const view_coords &c) const override {
    return T::valid_for(c);
  }

private:
  const char *n;
  const char *f;
};
} // namespace detail

/*
  Creates a fractal which is really a sequence of different fractals which
  can be run at different levels of precision.
*/
template <typename... Ts>
detail::MultiPrecisionFactory<Ts...> make_fractal(const char *name,
                                                  const char *family) {
  return {name, family};
}

template <typename... Ts>
detail::MultiPrecisionFactory<Ts...> make_fractal(const char *name) {
  return {name, name};
}

}; // namespace fractals

