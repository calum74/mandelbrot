#include "fractal.hpp"
#include "mandelbrot.hpp"
#include "orbit.hpp"

class SimpleMandeldrop : public fractals::PointwiseCalculation {
public:
  using Real = double;
  using Complex = std::complex<Real>;

  static bool valid_for(const view_coords &c) {
    return c.r <= 3 && convert<double>(c.r) > 1e-2; // 10;
  }

  static view_coords initial_coords() { return {0, -1, 3, 500}; }

  SimpleMandeldrop(const view_coords &c, int w, int h, std::atomic<bool> &stop)
      : max_iterations(c.max_iterations), coords(c, w, h) {}

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
  const int max_iterations;
  const fractals::plane<Real> coords;
};

const fractals::PointwiseFractal &naiveMandeldrop =
    fractals::make_fractal<SimpleMandeldrop>("Mandeldrop (low precision)");

template <typename LowPrecisionComplex, typename HighPrecisionComplex,
          typename Calculation>
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

  using SmallReal = typename LowPrecisionComplex::value_type;
  using BigReal = typename HighPrecisionComplex::value_type;

  static HighPrecisionComplex map(const HighPrecisionComplex &z) {
    auto n = inverse(mandelbrot::norm(z));
    return {-imag(z) * n, -real(z) * n};
  }

  LowPrecisionComplex map_delta(const LowPrecisionComplex &d) const {
    // Note that we can calculate deltas using low precision (`double`) complex
    // numbers. We just need to rearrange our calculation of the delta to avoid
    // loss of precision.
    return LowPrecisionComplex{0, 1} * d / (c0 * (c0 + d));
  }

  // Initialize the fractal. We calculate the high precision reference orbit
  // at the center of the view. Because this calculation can be time-consuming,
  // we provide a "stop" flag which is used to exit the calculation early if the
  // view changes, to keep the UI responsive.
  PerturbatedMandeldropCalculation(const view_coords &c, int w, int h,
                                   std::atomic<bool> &stop)
      : max_iterations(c.max_iterations), coords(c, w, h), ref_x(w / 2),
        ref_y(h / 2), c0{convert<SmallReal>(coords.x0) + (coords.dx * ref_x),
                         convert<SmallReal>(coords.y0) + (coords.dy * ref_y)},
        reference_orbit{
            mandelbrot::make_basic_orbit<Calculation>(map(HighPrecisionComplex{
                coords.x0 + fractals::convert<BigReal>(coords.dx * ref_x),
                coords.y0 + fractals::convert<BigReal>(coords.dy * ref_y)})),
            max_iterations, stop} {}

  // Are the given coordinates valid. Use this to prevent zooming out too far
  // or to select a different implementation for different resolutions.
  // The call to `valid_precision` checks the size of the radius relative to
  // the size of a BigReal to make sure we have sufficient accuracy.
  static bool valid_for(const view_coords &c) {
    return c.r < 2 && valid_precision(BigReal{c.r});
  }

  // The initial coordinates to view the Mandelbrot set.
  // Also specify the number of iterations (500).
  static view_coords initial_coords() { return {0, -1, 3, 500}; }

  // Calculates a single point of the fractal, at position (x,y).
  // Look up the actual coordinates (or in this case, the delta from the center
  // (ref_x, ref_y)) from the plane.
  double calculate(int x, int y) const override {
    LowPrecisionComplex delta =
        map_delta({coords.dx * (x - ref_x), coords.dy * (y - ref_y)});

    // The function `make_relative_orbit` will skip some iterations,
    // use `z.iteration()` to find out which iteration we are on.
    auto z = reference_orbit.make_relative_orbit(delta, this->max_iterations);

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
  const int max_iterations;

  // A mapping from points in the image to points in the complex plane.
  const fractals::plane<typename HighPrecisionComplex::value_type,
                        typename LowPrecisionComplex::value_type>
      coords;

  // Where in the image the reference orbit is.
  // Currently always at the center of the image.
  const int ref_x, ref_y;

  // The position of the reference orbit.
  const LowPrecisionComplex c0;

  // The calculated reference orbit, together with Taylor series terms for the
  // epsilon/dz for each iteration.
  mandelbrot::stored_taylor_series_orbit<
      LowPrecisionComplex,
      mandelbrot::basic_orbit<HighPrecisionComplex, Calculation>>
      reference_orbit;
};

template <int N>
using MD = PerturbatedMandeldropCalculation<
    std::complex<double>, std::complex<fractals::high_precision_real<N>>,
    mandelbrot::mandelbrot_calculation<2>>;

const fractals::PointwiseFractal &mandeldrop_fractal =
    fractals::make_fractal<SimpleMandeldrop, MD<4>, MD<6>, MD<10>,
                           MD<16> /*, MB<20> */>("Mandeldrop");

// Not used
class SimpleCubicMandelbrot : public fractals::PointwiseCalculation {
public:
  using Real = double;
  using Complex = std::complex<Real>;

  static bool valid_for(const view_coords &c) { return c.r < 2; }

  static view_coords initial_coords() { return {0, 0, 2, 500}; }

  SimpleCubicMandelbrot(const view_coords &c, int w, int h,
                        std::atomic<bool> &stop)
      : max_iterations(c.max_iterations), coords(c, w, h) {}

  double calculate(int x, int y) const override {
    auto p = coords(x, y);
    Complex c{p.x, p.y};
    Complex z = 0;
    int i = 0;
    while (!mandelbrot::escaped(z)) {
      if (i++ >= max_iterations)
        return 0;
      z = z * z * z + c;
    }

    return i;
  }

private:
  const int max_iterations;
  const fractals::plane<Real> coords;
};

const fractals::PointwiseFractal &simple_cubic_mandelbrot_fractal =
    fractals::make_fractal<SimpleCubicMandelbrot>(
        "Cubic Mandelbrot (low precision)");