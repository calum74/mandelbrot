#pragma once
#include "pointwise_fractal.hpp"
#include "view_coords.hpp"

#include <atomic>

namespace fractals {

/*
  The calculation of a fractal, where we can compute each (x,y) point
  independently and perhaps in parallel.

  A virtual function call per point may seem excessive, but it's usually
  negligible overhead relative to computing each point.
*/
class pointwise_calculation {
public:
  using view_coords = fractals::view_coords;
  virtual ~pointwise_calculation() = default;

  // Compute the point at the pixel position (x,y).
  virtual double calculate(int x, int y) const = 0;

  // Gets the average number of iterations, for fractals where this makes sense.
  virtual double average_iterations() const;

  // Gets the average number of skipped iterations, for fractals where this
  // makes sense.
  virtual double average_skipped() const;

  virtual void initialize(const view_coords &c, int x, int y,
                          std::atomic<bool> &stop) = 0;
};

/*
  A fractal that can be calculated in a point-wise fashion.
*/
class pointwise_calculation_factory {
public:
  virtual ~pointwise_calculation_factory() = default;

  virtual std::string name() const = 0;

  virtual std::string family() const = 0;

  // Creates and initializes a new calculation.
  // Called once for each generated fractal image.
  virtual std::shared_ptr<pointwise_calculation>
  create(const view_coords &c, int x, int y, std::atomic<bool> &stop) const = 0;

  // Retrieve the initial coordinates for this fractal.
  virtual view_coords initial_coords() const = 0;

  // Queries whether the given coordinates are valid for this fractal type.
  virtual bool valid_for(const view_coords &c) const = 0;
};

} // namespace fractals