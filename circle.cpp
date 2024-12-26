#include "fractal.hpp"

// Define a new fractal by inheriting from fractals::PointwiseCalculation
// In this case we'll specify the data-type of real numbers as a template
// parameter.
template <typename Real> class Circle : public fractals::PointwiseCalculation {
public:
  // Your coordinate system for this fractal.
  // Real specifies how you're going to store your coordinates.
  fractals::plane<Real> coords;

  // The constructor sets the viewing plane.
  // ViewCoords contains high-precision coordinates, so you'll want to translate
  // this into your local coordinate system. We can do some initialization at
  // this point, such as calculating reference orbits. Initialization can be
  // time-consuming, so we'll provide a `stop` cancellation token. The
  // constructor should regularly check `stop` and exit early to keep the UI
  // responsive.
  Circle(const ViewCoords &c, int w, int h, std::atomic<bool> &stop)
      : coords(c, w, h) {}

  // Returns whether the given fractal
  static bool valid_for(const ViewCoords &c) { return true; }

  static ViewCoords initial_coords() { return {0, 0, 1.5, 2}; }

  double calculate(int x, int y) const override {
    auto point = coords(x, y);

    return point.x * point.x + point.y * point.y >= 1.0;
  }
};

const fractals::PointwiseFractal &circle_fractal =
    fractals::make_fractal<Circle<double>,
                           Circle<fractals::high_precision_real<6>>>("Circle");
