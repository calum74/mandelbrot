#pragma once
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
this each iteration. size: The dimensions (width,height) of this region as a
complex number. central_orbit: A stored orbit generated from high-precision
numbers. starts at {0,0} stop: A cancellation token

    iteration is initially 0.
    j is initially 0.

*/
template <Complex LowPrecisionType, Complex DeltaType, Complex TermType,
          unsigned long Terms, Calculation Calc,
          RandomAccessOrbit HighPrecisionReferenceOrbit, typename OutputFn>
void magic(int iteration, int max_iterations, int x0, int y0, int w, int h,
           DeltaType epsilon, DeltaType size,
           const HighPrecisionReferenceOrbit &central_orbit, int j,
           std::atomic<bool> &stop, OutputFn fn) {

  auto size_norm = fractals::norm(size); // !! Should be previous size_norm/4

  // Start a new Taylor series for the epsilon relative to the current orbit
  // (not the central orbit)
  std::array<TermType, Terms> epsilon_terms;

  LowPrecisionType z;

  do {
    z = central_orbit[j];

    epsilon_terms = Calc::step_delta(z + epsilon, epsilon_terms);
    epsilon = Calc::step_epsilon(z, epsilon);
    epsilon = fractals::normalize(epsilon);
    j++;

    // TODO: Reference magic here

    iteration++;

    if (stop)
      return;
  } while (!stop && iteration < max_iterations && !escaped(z) &&
           size_norm < maximum_norm(epsilon_terms));

  if (iteration < max_iterations) {

    // Branch step - split the current orbit into 4 parts
    auto half_size = size * 0.5;
    auto quarter_size = half_size * 0.5;

    // For the top-left quadrant
    // Calculate the terms

    auto upper_left_center = -quarter_size;

    auto upper_left_epsilon =
        epsilon + evaluate_epsilon(upper_left_center, epsilon_terms);

    magic(iteration, max_iterations, x0, y0, w / 2, h / 2, upper_left_epsilon,
          half_size, central_orbit, stop);
  }
}
} // namespace mandelbrot