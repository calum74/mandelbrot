// This is an implementation of the classic Mandelbrot set fractal.

#include "mandelbrot.hpp"
#include "fractal.hpp"
#include "high_exponent_real.hpp"
#include "orbit.hpp"
#include "orbit_manager.hpp"

// Calculate the Mandelbrot set using perturbations and
// Taylor series to skip iterations. The algorithms are implemented in
// orbit.hpp. The class is templated so that we can configure the data types
// for higher resolution rendering if required.
template <mandelbrot::Complex LowPrecisionType, mandelbrot::Complex DeltaType,
          mandelbrot::Complex TermType, mandelbrot::Complex HighPrecisionType,
          mandelbrot::Calculation Calculation, int Terms, int TermPrecision>
class PerturbatedMandelbrotCalculation : public fractals::PointwiseCalculation {
public:
  using SmallReal = typename LowPrecisionType::value_type;
  using DeltaReal = typename DeltaType::value_type;
  using HighPrecisionReal = typename HighPrecisionType::value_type;

  using reference_orbit_type =
      mandelbrot::basic_orbit<HighPrecisionType, Calculation>;

  HighPrecisionType center;

  void initialize(const view_coords &c, int w, int h,
                  std::atomic<bool> &stop) override {

    max_iterations = c.max_iterations;

    auto new_center = HighPrecisionType{c.x, c.y};

    auto delta = convert<DeltaType>(new_center - center);

    center = new_center;

    fractals::plane<HighPrecisionReal, DeltaReal> new_coords(c, w, h);

    coords = new_coords;
    ref_x = w / 2;
    ref_y = h / 2;

    if (stop)
      return;

    reference_orbit_type init(HighPrecisionType{
        coords.x0 + fractals::convert<HighPrecisionReal>(coords.dx) *
                        fractals::convert<HighPrecisionReal>(ref_x),
        coords.y0 + fractals::convert<HighPrecisionReal>(coords.dy) *
                        fractals::convert<HighPrecisionReal>(ref_y)});

    orbits.new_view(
        delta, DeltaType{DeltaReal(0.5) * coords.w, DeltaReal(0.5) * coords.h},
        3, init, max_iterations, stop);
  }

  // Initialize the fractal. We calculate the high precision reference orbit
  // at the center of the view. Because this calculation can be time-consuming,
  // we provide a "stop" flag which is used to exit the calculation early if the
  // view changes, to keep the UI responsive.
  PerturbatedMandelbrotCalculation(const view_coords &c, int w, int h,
                                   std::atomic<bool> &stop) {}

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
    DeltaType delta = {coords.dx * DeltaReal(x - ref_x),
                       coords.dy * DeltaReal(y - ref_y)};

    auto z = orbits.lookup(delta, max_iterations);

    points_calculated++;
    skipped_iterations += z.iteration();
    while (fractals::norm(*z) <= SmallReal(1 << 16)) {
      if (z.iteration() >= this->max_iterations)
        return 0;
      ++z;
    }

    total_iterations += z.iteration();

    // This calculation creates a "fractional" iteration
    // used for smoother rendering.
    auto zn = log(fractals::norm(*z)) / 2;
    auto nu = log(log(zn) / log(2)) / log(2);
    return z.iteration() + 1 - nu;
  }

private:
  // The maxumum number of iterations / bailout value.
  int max_iterations;

  // A mapping from points in the image to points in the complex plane.
  fractals::plane<HighPrecisionReal, DeltaReal> coords;

  // Where in the image the reference orbit is.
  // Currently always at the center of the image.
  int ref_x, ref_y;

  // The calculated reference orbit, together with Taylor series terms for the
  // epsilon/dz for each iteration.
  mandelbrot::orbit_manager<LowPrecisionType, DeltaType, TermType, Terms,
                            TermPrecision, reference_orbit_type>
      orbits;
};

double fractals::PointwiseCalculation::average_iterations() const { return 0; }

double fractals::PointwiseCalculation::average_skipped() const { return 0; }

void fractals::PointwiseCalculation::initialize(const view_coords &c, int x,
                                                int y,
                                                std::atomic<bool> &stop) {}

template <int N, int P, int T = 4, int Tolerance = 100,
          typename DeltaType = std::complex<double>>
using MB = PerturbatedMandelbrotCalculation<
    std::complex<double>, DeltaType,
    std::complex<fractals::high_exponent_real<double>>,
    std::complex<fractals::high_precision_real<P>>,
    mandelbrot::mandelbrot_calculation<N>, T, Tolerance>;

template <int N, int P, int T = 4, int Tolerance = 100>
using MB_high =
    MB<N, P, T, Tolerance, std::complex<fractals::high_exponent_real<double>>>;

// TODO: a more elegant way to choose parameters

// Supply a list of fractals to `make_fractal`, which will create a factory
// that selects the best fractal at each resolution. We need different
// implementations at different resolutions so that we don't lose precision or
// use a slower algorithm than necessary.
// For the highest precition, we need to increase the precision on the epsilon
// as well
const fractals::PointwiseFractal &mandelbrot_fractal =
    fractals::make_fractal<MB<2, 3, 4, 1000>, MB<2, 6>, MB<2, 10>, MB<2, 18>,
                           MB_high<2, 32>, MB_high<2, 40>, MB_high<2, 64>>(
        "Mandelbrot (power 2)", "mandelbrot");

// Cubic Mandelbrot has no glitches with 3 Taylor series terms, but
// glitches quite badly with 4 terms. On the other hand, Square mandelbrot works
// better with 4 terms.
const fractals::PointwiseFractal &mandelbrot3_fractal =
    fractals::make_fractal<MB<3, 4, 3, 10000>, MB<3, 6, 3, 10000>,
                           MB<3, 10, 3, 10000>, MB<3, 18, 3, 10000>>(
        "Cubic Mandelbrot (power 3)");

const fractals::PointwiseFractal &mandelbrot4_fractal =
    fractals::make_fractal<MB<4, 4>, MB<4, 6>, MB<4, 10>, MB<4, 18>>(
        "Mandelbrot (power 4)");

const fractals::PointwiseFractal &mandelbrot5_fractal =
    fractals::make_fractal<MB<5, 4>, MB<5, 6>, MB<5, 10>, MB<5, 18>>(
        "Mandelbox (power 5)");

const fractals::PointwiseFractal &mandelbrot6_fractal =
    fractals::make_fractal<MB<6, 4>, MB<6, 6>, MB<6, 10>, MB<6, 18>>(
        "Mandelbrot (power 6)");

const fractals::PointwiseFractal &mandelbrot7_fractal =
    fractals::make_fractal<MB<7, 4>, MB<7, 6>, MB<7, 10>, MB<7, 18>>(
        "Mandelflake (power 7)");
