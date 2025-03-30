#pragma once
#include "fractal.hpp"
#include "view_coords.hpp"

#include <atomic>
#include <vector>

namespace fractals {

struct orbital_point {
  int x, y;
  int iteration;
};

using displayed_orbit = std::vector<orbital_point>;

/*
  The calculation of a fractal, where we can compute each (x,y) point
  independently and perhaps in parallel.

  A virtual function call per point may seem excessive, but it's usually
  negligible overhead relative to computing each point.
*/
class fractal_calculation {
public:
  using view_coords = fractals::view_coords;
  virtual ~fractal_calculation() = default;

  // Compute the point at the pixel position (x,y).
  virtual double calculate(int x, int y) const = 0;

  // Gets the average number of iterations, for fractals where this makes sense.
  virtual double average_iterations() const;

  // Gets the average number of skipped iterations, for fractals where this
  // makes sense.
  virtual double average_skipped() const;

  virtual void initialize(const view_coords &c, int x, int y,
                          std::atomic<bool> &stop) = 0;

  // Gets the pixels in an individual orbit.
  // The output will be sparse because it will only output points
  // that are within the current view.
  virtual void get_orbit(int x, int y,
                         std::vector<orbital_point> &points_out, std::atomic<bool> & stop) const;
};

/*
  A fractal that can be calculated in a point-wise fashion.
*/
class fractal_calculation_factory {
public:
  virtual ~fractal_calculation_factory() = default;

  virtual std::string name() const = 0;

  virtual std::string family() const = 0;

  // Creates and initializes a new calculation.
  // Called once for each generated fractal image.
  virtual std::shared_ptr<fractal_calculation>
  create(const view_coords &c, int x, int y, std::atomic<bool> &stop) const = 0;

  // Retrieve the initial coordinates for this fractal.
  virtual view_coords initial_coords() const = 0;

  // Queries whether the given coordinates are valid for this fractal type.
  virtual bool valid_for(const view_coords &c) const = 0;
};

} // namespace fractals