/*
  Sample showing how to add new fractals.
*/

#include "fractal.hpp"

// Define a new fractal by inheriting from fractals::PointwiseCalculation
// In this case we'll specify the data-type of real numbers as a template
// parameter.
template <typename Real> class Circle : public fractals::PointwiseCalculation {
public:
  // Your coordinate system for this fractal.
  // Real specifies how you're going to store your coordinates.
  fractals::plane<Real> coords;

  void initialize(const view_coords &c, int w, int h,
                  std::atomic<bool> &stop) override {
    coords = {c, w, h};
  }

  // Returns whether the given fractal
  static bool valid_for(const view_coords &c) { return true; }

  static view_coords initial_coords() { return {0, 0, 1.5, 2}; }

  double calculate(int x, int y) const override {
    auto point = coords(x, y);

    return point.x * point.x + point.y * point.y >= 1.0;
  }
};

const fractals::PointwiseCalculationFactory &circle_fractal =
    fractals::make_fractal<Circle<double>,
                           Circle<fractals::high_precision_real<6>>>("Circle");
