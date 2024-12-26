#pragma once

#include <memory>

#include "plane.hpp"
#include "view_coords.hpp"

/* Basic infrastructure for rendering fractals. */

namespace fractals {

class PointwiseCalculation {
public:
  using view_coords = fractals::view_coords;
  virtual ~PointwiseCalculation() = default;

  virtual double calculate(int x, int y) const = 0;

  virtual double average_iterations() const;
  virtual double average_skipped() const;
};

class PointwiseFractal {
public:
  virtual const char *name() const = 0;
  virtual std::unique_ptr<PointwiseCalculation>
  create(const view_coords &c, int x, int y, std::atomic<bool> &stop) const = 0;
  virtual view_coords initial_coords() const = 0;
  virtual bool valid_for(const view_coords &c) const = 0;
};

namespace detail {

template <typename... Ts> class MultiPrecisionFactory;

template <typename T, typename... Ts>
class MultiPrecisionFactory<T, Ts...> : public PointwiseFractal {
public:
  MultiPrecisionFactory(const char *name) : n{name}, tail{name} {}

  view_coords initial_coords() const override { return T::initial_coords(); }
  const char *name() const override { return n; }

  std::unique_ptr<PointwiseCalculation>
  create(const view_coords &c, int x, int y,
         std::atomic<bool> &stop) const override {
    return T::valid_for(c) ? std::make_unique<T>(c, x, y, stop)
                           : tail.create(c, x, y, stop);
  }

  bool valid_for(const view_coords &c) const override {
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
  create(const view_coords &c, int x, int y,
         std::atomic<bool> &stop) const override {
    return std::make_unique<T>(c, x, y, stop);
  }

  view_coords initial_coords() const override { return T::initial_coords(); }

  const char *name() const override { return n; }

  bool valid_for(const view_coords &c) const override {
    return T::valid_for(c);
  }

private:
  const char *n;
};
} // namespace detail

template <typename... Ts>
detail::MultiPrecisionFactory<Ts...> make_fractal(const char *name) {
  return {name};
}

}; // namespace fractals

