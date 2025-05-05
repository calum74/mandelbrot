// This is an implementation of the classic Mandelbrot set fractal.

#include "mandelbrot.hpp"
#include "generate_fractal.hpp"
#include "orbit.hpp"
#include "orbit_manager.hpp"

// Calculate the Mandelbrot set using perturbations and
// Taylor series to skip iterations. The algorithms are implemented in
// orbit.hpp. The class is templated so that we can configure the data types
// for higher resolution rendering if required.
template <numbers::complex LowPrecisionType, numbers::complex DeltaType,
          numbers::complex TermType, numbers::complex HighPrecisionType,
          mandelbrot::Calculation Calculation, int Terms, int TermPrecision1,
          int TermPrecision2, int NumOrbits>
class PerturbatedMandelbrotCalculation : public fractals::fractal_calculation {
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

    auto new_center = HighPrecisionType{number_cast<HighPrecisionReal>(c.x),
                                        number_cast<HighPrecisionReal>(c.y)};

    auto delta = numbers::number_cast<DeltaType>(new_center - center);

    center = new_center;

    fractals::plane<HighPrecisionReal, DeltaReal> new_coords(c, w, h);

    coords = new_coords;
    ref_x = w / 2;
    ref_y = h / 2;

    if (stop)
      return;

    reference_orbit_type init(HighPrecisionType{
        coords.x0 + numbers::number_cast<HighPrecisionReal>(coords.dx) *
                        numbers::number_cast<HighPrecisionReal>(ref_x),
        coords.y0 + numbers::number_cast<HighPrecisionReal>(coords.dy) *
                        numbers::number_cast<HighPrecisionReal>(ref_y)});

    orbits.new_view(
        delta, DeltaType{DeltaReal(0.5) * coords.w, DeltaReal(0.5) * coords.h},
        NumOrbits, init, max_iterations, stop);
  }

  // Are the given coordinates valid. Use this to prevent zooming out too far
  // or to select a different implementation for different resolutions.
  // The call to `valid_precision` checks the size of the radius relative to
  // the size of a BigReal to make sure we have sufficient accuracy.
  static bool valid_for(const view_coords &c) {
    return c.r <= 2 &&
           numbers::valid_precision(number_cast<HighPrecisionReal>(c.r));
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

    auto z = orbits.lookup(delta, max_iterations, true);

    points_calculated++;
    skipped_iterations += z.iteration();
    while (numbers::norm(*z) <= SmallReal(1 << 16)) {
      if (z.iteration() >= this->max_iterations)
        return 0;
      ++z;
    }

    total_iterations += z.iteration();

    // This calculation creates a "fractional" iteration
    // used for smoother rendering.
    auto zn = log2(numbers::norm(*z));
    auto nu = log2(zn) / log2(Calculation::order);
    return z.iteration() + 1 - nu;
  }

  void get_orbit(int x, int y, fractals::displayed_orbit &output,
                 std::atomic<bool> &stop) const override {

    mandelbrot::basic_orbit<HighPrecisionType, Calculation> orbit(
        HighPrecisionType{coords.get_x(x), coords.get_y(y)});

    output.clear();
    for (int iteration = 0;
         !mandelbrot::escaped(*orbit) && iteration < max_iterations && !stop;
         ++iteration) {
      ++orbit;

      auto delta = numbers::number_cast<DeltaType>(orbit.distance());
      auto max_norm = coords.w * coords.h; // Not correct

      if (numbers::norm(delta) < max_norm) {
        auto p = *orbit;
        auto x = coords.to_x(numbers::real_part(p));
        auto y = coords.to_y(numbers::imag_part(p));
        if(x>=0 && y>=0)
          output.push_back({x,y, iteration});
      }
    }
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
                            TermPrecision1, TermPrecision2,
                            reference_orbit_type>
      orbits;
};

// Configure a Mandelbrot set fractal based on its power and precision.
template <int Power> struct mandelbrot_generator {
  template <int Precision> struct precision {

    // The different numerical types we need.
    // Uses the helper type `fractals::complex_number` to specify a complex
    // number of the desired precision.
    using low_precision_type = numbers::complex_number<48, 0, 0>;
    using delta_type = numbers::complex_number<0, -2 * Precision, 0>;
    using term_type = numbers::complex_number<0, -1000000, 1000000>;
    using high_precision_type = numbers::complex_number<Precision, 0, 0>;

    // The equations we need
    using calculation_type = mandelbrot::mandelbrot_calculation<Power>;

    // Various fudge factors. The reason these are here is to render
    // more quickly without introducing glitches, otherwise we would
    // just calculate everything to the highest precision but that would be
    // too slow.
    static constexpr int series_terms = 4;
    static constexpr int term_precision1 = 25;
    static constexpr int term_precision2 = 100;
    static constexpr int num_orbits = Precision <= 1024 ? 3 : 1;

    using type = PerturbatedMandelbrotCalculation<
        low_precision_type, delta_type, term_type, high_precision_type,
        calculation_type, series_terms, term_precision1, term_precision2,
        num_orbits>;
  };
};

const fractals::fractal &mandelbrot_fractal =
    fractals::generate_fractal<mandelbrot_generator<2>>("Mandelbrot set",
                                                        "mandelbrot");

const fractals::fractal &mandelbrot3_fractal =
    fractals::generate_fractal<mandelbrot_generator<3>>(
        "Cubic Mandelbrot (power 3)");

const fractals::fractal &mandelbrot4_fractal =
    fractals::generate_fractal<mandelbrot_generator<4>>("Mandelbrot (power 4)");

const fractals::fractal &mandelbrot5_fractal =
    fractals::generate_fractal<mandelbrot_generator<5>>("Mandelbox (power 5)");

const fractals::fractal &mandelbrot6_fractal =
    fractals::generate_fractal<mandelbrot_generator<6>>("Mandelbrot (power 6)");

const fractals::fractal &mandelbrot7_fractal =
    fractals::generate_fractal<mandelbrot_generator<7>>("Mandelbrot (power 7)");
