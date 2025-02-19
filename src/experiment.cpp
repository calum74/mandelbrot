
#include "bla.hpp"
#include "high_exponent_real.hpp"
#include "mandelbrot.hpp"
#include "make_fractal.hpp"

#include <mutex>
#include <optional>

// Calculate the Mandelbrot set using perturbations and
// Taylor series to skip iterations. The algorithms are implemented in
// orbit.hpp. The class is templated so that we can configure the data types
// for higher resolution rendering if required.
template <mandelbrot::Complex LowPrecisionType, mandelbrot::Complex DeltaType,
          mandelbrot::Complex TermType, mandelbrot::Complex HighPrecisionType,
          mandelbrot::Calculation Calculation, int Terms, int TermPrecision>
class ExperimentalMandelbrotCalculation
    : public fractals::pointwise_calculation {
public:
  using SmallReal = typename LowPrecisionType::value_type;
  using DeltaReal = typename DeltaType::value_type;
  using HighPrecisionReal = typename HighPrecisionType::value_type;

  using reference_type =
      mandelbrot::basic_orbit<HighPrecisionType, Calculation>;
  using stored_type =
      mandelbrot::stored_orbit<LowPrecisionType, reference_type>;

  void initialize(const view_coords &c, int w, int h,
                  std::atomic<bool> &stop) override {
    coords = {c, w, h};
    max_iterations = c.max_iterations;

    reference_type reference_orbit(HighPrecisionType{c.x, c.y});
    std::lock_guard<std::mutex> lock(m);
    stored_orbit = {reference_orbit, c.max_iterations, stop};
    orbit = {*stored_orbit};
    skipped_iterations = 0;
    points_calculated = 0;
    total_iterations = 0;
  }

  std::optional<stored_type> stored_orbit;
  // mutable mandelbrot::linear_orbit<DeltaType, stored_type> orbit;
  mutable mandelbrot::bilinear_orbit<DeltaType, TermType, stored_type> orbit;

  // Are the given coordinates valid. Use this to prevent zooming out too far
  // or to select a different implementation for different resolutions.
  // The call to `valid_precision` checks the size of the radius relative to
  // the size of a BigReal to make sure we have sufficient accuracy.
  static bool valid_for(const view_coords &c) {
    // !! Bug in comparison here - allows up to 3.
    return c.r <= 2 && valid_precision(HighPrecisionReal{c.r});
  }

  // The initial coordinates to view the Mandelbrot set.
  // Also specify the number of iterations (500).
  static view_coords initial_coords() {
    return {Calculation::order > 2 ? 0.0 : -0.5, 0.0, 2, 500};
  }

  mutable long skipped_iterations = 0;
  mutable long points_calculated = 0;
  mutable long total_iterations = 0;

  double average_iterations() const override {
    return double(total_iterations) / points_calculated;
  }

  double average_skipped() const override {
    return double(skipped_iterations) / points_calculated;
  }

  mutable std::atomic<int> iterations_skipped;

  mutable std::mutex m;
  // Calculates a single point of the fractal, at position (x,y).
  // Look up the actual coordinates (or in this case, the delta from the
  // center (ref_x, ref_y)) from the plane.
  double calculate(int x, int y) const override {

    DeltaType delta = {coords.dx * (x - coords.pw / 2),
                       coords.dy * (y - coords.ph / 2)};

    std::lock_guard<std::mutex> lock(m);
    auto r = orbit.get(delta, max_iterations);
    ++points_calculated;
    skipped_iterations += orbit.get_skipped_iterations();
    total_iterations += r;
    return r;
  }

private:
  // A mapping from points in the image to points in the complex plane.
  fractals::plane<HighPrecisionReal, DeltaReal> coords;
  int max_iterations;
};

template <int N, int P, int T = 4, int Tolerance = 100,
          typename DeltaType = std::complex<double>>
using MBX = ExperimentalMandelbrotCalculation<
    std::complex<double>, DeltaType,
    std::complex<fractals::high_exponent_double>,
    std::complex<fractals::high_precision_real<P>>,
    mandelbrot::mandelbrot_calculation<N>, T, Tolerance>;

// Nothing here now
const fractals::pointwise_fractal &experimental_fractal =
    fractals::make_fractal<MBX<2, 3*64>, MBX<2, 6*64>, MBX<2, 10*64>, MBX<2, 18*64>>(
        "Experiment", "mandelbrot");
