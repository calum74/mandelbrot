/*
    Welcome to the tutorial for this library!

    This library contains many bits and pieces for implementing
    the Mandelbrot set.
*/

#include "mandelbrot.hpp"
#include "orbit.hpp"
#include <iomanip>

int main() {

  // To create a basic fractal, we'll need a fractal factory
  // For example, the fractal for the Mandelbrot set is called
  // `mandelbrot_fractal`, defined in `mandelbrot.hpp`

  // Get the initial coordinates for this fractal.
  // This also specifies the maximum number of iterations.
  auto coords = mandelbrot_fractal.initial_coords();

  // Create a "calculation" for the fractal, supplying the dimensions to
  // calculate and also a `stop` token (a cancellation token) in case we want to
  // cancel the calculation at any time.
  std::atomic<bool> stop;
  auto view = mandelbrot_fractal.create(coords, 100, 100, stop);

  // We can then calculate any point we like.
  // It's up to the caller to decide on the evaluation order and
  // whether to calculate these in threads.
  auto p = view->calculate(25, 25);
  std::cout << "Escaped after " << p << " iterations\n";

  // We can also access algorithms

  // A "basic" orbit is an iterator that just evaluates each point z in the
  // sequence. Use * to get the current value, and ++ to advance to the next
  // point.

  // To create a basic orbit, we'll supply a data-type for each point in the
  // orbit (for example, std::complex<double>), and specify which calculation we
  // are using. For example, `mandelbrot_calculation<2>` computes the standard
  // Mandelbrot set.
  mandelbrot::basic_orbit<std::complex<double>,
                          mandelbrot::mandelbrot_calculation<2>>
      orbit1({0.5, 0.5});

  int iterations = 0;
  while (!mandelbrot::escaped(*orbit1)) {
    ++iterations;
    ++orbit1;
  }

  // Escapes after 5 iterations
  std::cout << "Escaped after " << iterations << " iterations\n";

  // To get the deep zooms, high precision arithmetic is needed.
  // A `double` datatype is only good up to about a radius of 1e-14.
  // Instead, we have implemented a very basic data-type for high-precision
  // arithmetic, where the template parameter specifies the number of 64-bit
  // integers used to store the number.
  using R = fractals::high_precision_real<8>;

  // You can convert high_precision_real to and from strings using iostreams
  std::stringstream ss{"1.234567890123456789012345"};
  R r1;
  ss >> r1;
  std::cout << std::setprecision(20) << r1 << std::endl;

  // You can supply high_precision_real<> numbers to std::complex,
  // and create basic (high-precision) orbits

  auto reference_orbit =
      mandelbrot::make_basic_orbit<mandelbrot::mandelbrot_calculation<2>>(
          std::complex<R>{0.49, 0.49});

  // In order to speed up Mandelbrot calculations, we can use
  // perturbation theory to perform most of our calculations using
  // low-precision arithmetic (like `std::complex<double>`), and only compute
  // one reference orbit.
  // The stored orbit only needs to store low-precision numbers like
  // std::complex<double>

  // Let's create a stored orbit, which lazily stores all of the complex numbers
  // in the orbit. Although the numbers stored in the orbit are of low precision
  // std::complex<double>, they have been derived from a high-precision orbit.
  auto stored_orbit = mandelbrot::make_stored_orbit<std::complex<double>>(
      reference_orbit, 100, stop);

  // Now, create a "relative orbit", which is a low-precision delta {0.01,0.01},
  // relative to our reference orbit.

  // Now we can iterate the relative orbit as if were a normal orbit.

  // The advanced algorithms can use Taylor series to skip large numbers of
  // iterations

  // Compute the Taylor series coefficients for the given reference orbit.
  mandelbrot::stored_taylor_series_orbit<
      std::complex<double>, std::complex<double>, std::complex<double>,
      mandelbrot::basic_orbit<std::complex<R>,
                              mandelbrot::mandelbrot_calculation<2>>,
      3, 100>
      taylor_series{reference_orbit, 100, stop};

  // Construct a relative orbit to the reference orbit.
  // This automatically skips some iterations in the orbit, given by
  // iteration() to give you the actually iteration number.
  auto relative =
      taylor_series.make_relative_orbit({0.01, 0.01}, 100, iterations);
  iterations = relative.iteration();

  std::cout << "Skipped " << iterations << " iterations\n";

  // Iterate the relative orbit as before, using perturbation.
  while (!mandelbrot::escaped(*relative) && iterations < 100) {
    ++iterations;
    ++relative;
  }

  std::cout << "Escaped after " << iterations << " iterations\n";

  // Let's look at Taylor series "clusters"

  // To create a cluster, we need to configure the data-types
  // orbit_type
  // delta_type
  // term_type
  // reference orbit
  // number of taylor-series terms (3)
  // Taylor-series precision factor.

  mandelbrot::taylor_series_cluster<
      std::complex<double>, std::complex<double>, std::complex<double>,
      mandelbrot::basic_orbit<std::complex<R>,
                              mandelbrot::mandelbrot_calculation<2>>,
      3, 100>
      cluster1{reference_orbit, 100, stop};

  // The constructor has evaluated the reference orbit to the
  // specified number of iterations (or, until it has escaped, whichever is
  // sooner).

  std::cout << "Reference orbit contains " << cluster1.primary_orbit.size()
            << " points\n";

  // We can then create a secondary orbit relative within the cluster
  auto secondary = cluster1.make_secondary_orbit({0.01, 0.0}, 100, stop);

  std::cout << "Secondary orbit contains " << secondary.size() << " points\n";

  iterations = 0;
  // And then our final orbit
  auto tertiary = secondary.make_relative_orbit({0.0, 0.01}, 100, iterations);

  // We can then skip iterations on tertiary
  std::cout << "Skipped " << iterations << " iterations\n";

  // Iterate the tertiary orbit as before, using perturbation.
  while (!mandelbrot::escaped(*tertiary) && iterations < 100) {
    ++iterations;
    ++tertiary;
  }

  std::cout << "Escaped after " << iterations << " iterations\n";
}