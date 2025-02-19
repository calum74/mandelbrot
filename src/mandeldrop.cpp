#include "fractal.hpp"
#include "mandelbrot.hpp"
#include "mandelbrot_adaptor.hpp"
#include "orbit.hpp"
#include "convert.hpp"

class SimpleMandeldrop : public fractals::PointwiseCalculation {
public:
  using Real = double;
  using Complex = std::complex<Real>;

  static bool valid_for(const view_coords &c) {
    return c.r <= 4 && convert<double>(c.r) > 1e-2; // 10;
  }

  static view_coords initial_coords() { return {0, -1, 3, 500}; }

  void initialize(const view_coords &c, int w, int h,
                  std::atomic<bool> &stop) override {
    max_iterations = c.max_iterations;
    coords = {c, w, h};
  }

  double calculate(int x, int y) const override {
    auto p = coords(x, y);
    Complex c0{p.x, p.y};

    // A conformal mapping, so this just projects the Mandelbrot set to a
    // different location, and locally it still looks like a Mandelbrot set.
    Complex c = Complex{0, -1} / c0;
    Complex z = 0;
    int i = 0;
    while (!mandelbrot::escaped(z)) {
      if (i++ >= max_iterations)
        return 0;
      z = z * z + c;
    }

    return i;
  }

private:
  int max_iterations;
  fractals::plane<Real> coords;
};

const fractals::PointwiseFractal &naiveMandeldrop =
    fractals::make_fractal<SimpleMandeldrop>("Mandeldrop (low precision)");

template <typename Adaptor, mandelbrot::Complex LowPrecisionType,
          mandelbrot::Complex DeltaType, mandelbrot::Complex TermType,
          mandelbrot::Complex HighPrecisionType,
          mandelbrot::Calculation Calculation, int Terms, int TermPrecision1,
          int TermPrecision2>
class PerturbatedMandeldropCalculation : public fractals::PointwiseCalculation {
public:
  /*
    The "Mandeldrop" is a transformation of the Mandelbrot set under

      c' = -i/c

    The `-i` on the top rotates the image 90 degrees, and `1/c` essentially
    inverts the image. Using perturbation, we have

      z' + d'

    as before. This should give d in terms of d' as

      d' = id/(c(c+d))

    So the basic algorithm is that we transform the central point by
    -i/c, and transform the delta by id/(c(c+d)).

    Then we can plug the c and d into our original high-performance Mandelbrot
    algorithms.

    Mandeldrop looks like the Mandelbrot set at higher magnification because
    this is a conformal mapping.
  */

  using SmallReal = typename LowPrecisionType::value_type;
  using DeltaReal = typename DeltaType::value_type;
  using HighPrecisionReal = typename HighPrecisionType::value_type;

  void initialize(const view_coords &c, int w, int h,
                  std::atomic<bool> &stop) override {
    max_iterations = c.max_iterations;
    coords = {c, w, h};
    ref_x = (w / 2);
    ref_y = (h / 2);
    c0 = {convert<SmallReal>(coords.x0) + (coords.dx * ref_x),
          convert<SmallReal>(coords.y0) + (coords.dy * ref_y)};
    reference_orbit = {
        mandelbrot::make_basic_orbit<Calculation>(Adaptor::map(
            HighPrecisionType{coords.x0 + fractals::convert<HighPrecisionReal>(
                                              coords.dx * ref_x),
                              coords.y0 + fractals::convert<HighPrecisionReal>(
                                              coords.dy * ref_y)})),
        max_iterations, stop};
  }

  // Are the given coordinates valid. Use this to prevent zooming out too far
  // or to select a different implementation for different resolutions.
  // The call to `valid_precision` checks the size of the radius relative to
  // the size of a BigReal to make sure we have sufficient accuracy.
  static bool valid_for(const view_coords &c) {
    return c.r <= 3 && valid_precision_for_inverse(HighPrecisionReal{c.r});
  }

  // The initial coordinates to view the Mandelbrot set.
  // Also specify the number of iterations (500).
  static view_coords initial_coords() { return {0, -1, 3, 500}; }

  mutable std::atomic<int> iterations_skipped;

  // Calculates a single point of the fractal, at position (x,y).
  // Look up the actual coordinates (or in this case, the delta from the center
  // (ref_x, ref_y)) from the plane.
  double calculate(int x, int y) const override {
    DeltaType delta = Adaptor::map_delta(
        DeltaType{coords.dx * (x - ref_x), coords.dy * (y - ref_y)}, c0);

    int skipped = iterations_skipped;
    // The function `make_relative_orbit` will skip some iterations,
    // use `z.iteration()` to find out which iteration we are on.
    auto z = reference_orbit.make_relative_orbit(delta, this->max_iterations,
                                                 skipped);
    iterations_skipped = skipped;

    while (norm(*z) <= (1 << 16)) {
      if (z.iteration() >= this->max_iterations)
        return 0;
      ++z;
    }

    // This calculation creates a "fractional" iteration
    // used for smoother rendering.
    auto zn = log(norm(*z)) / 2;
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

  // The position of the reference orbit.
  LowPrecisionType c0;

  // The calculated reference orbit, together with Taylor series terms for the
  // epsilon/dz for each iteration.
  mandelbrot::stored_taylor_series_orbit<
      LowPrecisionType, DeltaType, TermType,
      mandelbrot::basic_orbit<HighPrecisionType, Calculation>, Terms,
      TermPrecision1, TermPrecision2>
      reference_orbit;
};

template <int N, int P, int Terms = 4, int TP1 = 20, int TP2 = 100,
          typename DeltaType = std::complex<double>>
using MD = PerturbatedMandeldropCalculation<
    mandelbrot::mandeldrop_adaptor, std::complex<double>, DeltaType,
    std::complex<fractals::high_exponent_double>,
    std::complex<fractals::high_precision_real<P>>,
    mandelbrot::mandelbrot_calculation<2>, Terms, TP1, TP2>;

const fractals::PointwiseFractal &mandeldrop_fractal =
    fractals::make_fractal<MD<2, 4>, MD<2, 6>, MD<2, 10>,
                           MD<2, 16> /*, MB<20> */>("Mandeldrop");
