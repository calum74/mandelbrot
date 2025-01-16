
#include "fractal.hpp"
#include "high_exponent_real.hpp"
#include "magic_algorithm.hpp"
#include "mandelbrot.hpp"
#include "orbit.hpp"
#include "reference_orbit_manager.hpp"

// Calculate the Mandelbrot set using perturbations and
// Taylor series to skip iterations. The algorithms are implemented in
// orbit.hpp. The class is templated so that we can configure the data types
// for higher resolution rendering if required.
template <mandelbrot::Complex LowPrecisionType, mandelbrot::Complex DeltaType,
          mandelbrot::Complex TermType, mandelbrot::Complex HighPrecisionType,
          mandelbrot::Calculation Calculation, int Terms, int TermPrecision>
class ExperimentalMandelbrotCalculation
    : public fractals::PointwiseCalculation {
public:
  using SmallReal = typename LowPrecisionType::value_type;
  using DeltaReal = typename DeltaType::value_type;
  using HighPrecisionReal = typename HighPrecisionType::value_type;

  // Initialize the fractal. We calculate the high precision reference orbit
  // at the center of the view. Because this calculation can be time-consuming,
  // we provide a "stop" flag which is used to exit the calculation early if the
  // view changes, to keep the UI responsive.
  ExperimentalMandelbrotCalculation(const view_coords &c, int w, int h,
                                    std::atomic<bool> &stop)
      : coords(c, w, h) {
    auto diagonal_size = DeltaType(coords.dx * fractals::convert<DeltaReal>(w),
                                   coords.dy * fractals::convert<DeltaReal>(h));

    pw = w;
    experiment.resize(w * h);

    mandelbrot::basic_orbit<HighPrecisionType, Calculation> orbit1(
        HighPrecisionType{
            coords.x0 + fractals::convert<HighPrecisionReal>(coords.dx) *
                            fractals::convert<HighPrecisionReal>(w / 2),
            coords.y0 + fractals::convert<HighPrecisionReal>(coords.dy) *
                            fractals::convert<HighPrecisionReal>(h / 2)});
    mandelbrot::stored_orbit<
        LowPrecisionType,
        mandelbrot::basic_orbit<HighPrecisionType, Calculation>>
        stored(orbit1, c.max_iterations, stop);

    mandelbrot::magic<LowPrecisionType, DeltaType, TermType, 4>(
        c.max_iterations, 0, 0, w, h, diagonal_size, stored, stop,
        [&](int x, int y, double it, int reference_it) {
          // std::cout << "(" << x << "," << y << "," << it << ")";
          experiment.at(x + y * w) = it;
        });
  }

  int pw;
  std::vector<double> experiment;

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

  // Calculates a single point of the fractal, at position (x,y).
  // Look up the actual coordinates (or in this case, the delta from the
  // center (ref_x, ref_y)) from the plane.
  double calculate(int x, int y) const override {
    return experiment[x + y * pw];
  }

private:
  // A mapping from points in the image to points in the complex plane.
  const fractals::plane<HighPrecisionReal, DeltaReal> coords;
};

template <int N, int P, int T = 4, int Tolerance = 100,
          typename DeltaType = std::complex<double>>
using MBX = ExperimentalMandelbrotCalculation<
    std::complex<double>, DeltaType,
    std::complex<fractals::high_exponent_real<double>>,
    std::complex<fractals::high_precision_real<P>>,
    mandelbrot::mandelbrot_calculation<N>, T, Tolerance>;

// Nothing here now
const fractals::PointwiseFractal &experimental_fractal =
    fractals::make_fractal<MBX<2, 3, 4, 1000>, MBX<2, 6>, MBX<2, 10>,
                           MBX<2, 18>>("Experiment");
