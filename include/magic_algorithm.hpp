#pragma once
#include "convert.hpp"
#include "orbit.hpp"

namespace mandelbrot {

/*
This is an attempt to calculate the Mandelbrot set more efficiently. It works by
creating a recursive "tree" of reference orbits, and different reference orbits
can share the early terms in the reference orbit until it is determined that the
orbits need to branch.

We don't need to construct each tree explicitly, but rather follow the outline
of the tree using recursion.

Starting with a reference orbit, and a distance `epsilon` from that reference
orbit, we'll calculate a Taylor series for the relative epsilon around the point
`z+e`.

The "trick" is that we can use the divergence of the Taylor series to tell us
when to branch. We want to branch as late as possible as this allows us to share
more of the trees.

However, we don't need to store the terms of the Taylor series, we only need to
monitor if they diverge.

  Keep on iterating `iteration` until
   - We reach max_iterations
   - Our size is very small
   - The series escapes
   - The series begins to diverge

Arguments:

    iteration: The current iteration starting from 0. Increments.
    max_iterations: The bailout value, does not change.
    x0, y0: The position of the top-left pixel.
    w, h: The size of the region in pixels.
    epsilon: The position of this orbit relative to the central orbit. We update
        this each iteration.

    diagnonal_size: The dimensions (width,height) of this region as a
        complex number.

central_orbit: A stored orbit generated from high-precision
        numbers. starts at {0,0} stop: A cancellation token

    iteration is initially 0.
    j is initially 0.
*/
template <Complex LowPrecisionType, Complex DeltaType, Complex TermType,
          unsigned long Terms, RandomAccessOrbit HighPrecisionReferenceOrbit,
          typename OutputFn>
void magic(int max_iterations, int x0, int y0, int w, int h, DeltaType diagonal,
           const HighPrecisionReferenceOrbit &central_orbit,
           std::atomic<bool> &stop, OutputFn fn, DeltaType central_delta = {},
           int iteration = 0, int j = 0) {

  using calculation = typename HighPrecisionReferenceOrbit::calculation;
  // Novelty #1:
  // We use the norm of the delta to inform us how much we allow the series to
  // diverge
  typename TermType::value_type size_norm =
      fractals::norm(diagonal); // !! Should be previous size_norm/4

  // Novelty #2: Construct a Taylor series for the epsilon relative to the
  // current orbit (not the central orbit)
  std::array<TermType, Terms> epsilon_terms;

  DeltaType epsilon = {0, 0};

  LowPrecisionType z;

  do {
    z = central_orbit[j];

    epsilon_terms = calculation::delta_terms(
        z + convert<LowPrecisionType>(epsilon), epsilon_terms);
    epsilon = calculation::step_epsilon(z, epsilon, central_delta);
    epsilon = fractals::normalize(epsilon);
    j++;

    if (j == central_orbit.size() - 1 || escaped(central_orbit[j]) ||
        fractals::norm(z) < convert<typename LowPrecisionType::value_type>(
                                fractals::norm(epsilon))) {
      // We have exceeded the bounds of the current orbit
      // We need to reset the current orbit.
      // Thanks to
      // https://github.com/ImaginaFractal/Algorithms/blob/main/Perturbation.cpp
      // https://philthompson.me/2022/Perturbation-Theory-and-the-Mandelbrot-set.html
      epsilon = z;
      j = 0;
    }

    iteration++;

    if (stop)
      return;
  } while (!stop && iteration < max_iterations && !escaped(z) &&
           size_norm < maximum_delta_norm(epsilon_terms));

  std::cout << "Iteration to " << iteration << std::endl;

  if (iteration < max_iterations) {

    // Branch step - split the current orbit into 4 parts
    auto half_size = diagonal * DeltaType(0.5);
    auto quarter_size = half_size * DeltaType(0.5);

    auto lower_right_delta = quarter_size;
    auto upper_left_delta = -quarter_size;
    auto upper_right_delta =
        DeltaType(-upper_left_delta.real(), upper_left_delta.imag());
    auto lower_left_delta = -upper_right_delta;

    auto upper_left_epsilon =
        epsilon + evaluate_epsilon(upper_left_delta, epsilon_terms);
    auto lower_right_epsilon =
        epsilon + evaluate_epsilon(lower_right_delta, epsilon_terms);
    auto upper_right_epsilon =
        epsilon + evaluate_epsilon(upper_right_delta, epsilon_terms);
    auto lower_left_epsilon =
        epsilon + evaluate_epsilon(lower_left_delta, epsilon_terms);

    magic<LowPrecisionType, DeltaType, TermType, Terms>(
        max_iterations, x0, y0, w / 2, h / 2, half_size, central_orbit, stop,
        fn, central_delta + upper_left_delta, iteration, j);

    // etc
  }
}
} // namespace mandelbrot