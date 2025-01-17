#pragma once
#include "orbit.hpp"

namespace mandelbrot {

/*
A reference branch is similar to a reference orbit, however it is stored in
sections. Nearby branches can share a common root, and we only need to create
new branches when the current branch runs out of road.
*/
template <typename OrbitType, typename DeltaType, typename TermType,
          unsigned long Terms, typename HighPrecisionReferenceOrbit>
class reference_branch {

  struct series_entry {
    OrbitType z;
    std::array<TermType, Terms> terms;

    DeltaType epsilon(DeltaType delta) const {
      return evaluate_epsilon(delta, terms);
    }

    DeltaType evaluate_z(DeltaType delta) const { return z + epsilon(delta); }

    typename TermType::value_type maximum_delta_norm() const {
      // TODO: Tweak this when done
      return mandelbrot::maximum_delta_norm<10, 1000>(terms);
    }

    perturbation_orbit<OrbitType, DeltaType, HighPrecisionReferenceOrbit> orbit;
  };

public:
  // The iteration number
  // Iteration 0 is the very start of the sequence, and has value z(0) = 0.
  // Each branch covers a set number of iterations, and has the invariant
  // iteration = previous->iteration + previous->size()
  int iteration;

  // If there is a previous branch, this is it
  // Invariant: iteration==0 or previous_branch!=0
  std::shared_ptr<reference_branch> parent;

  // The delta of the orbit relative to the previous branch
  DeltaType delta_to_parent;

  reference_branch() : iteration(0), delta_to_parent(0) {}

  reference_branch(const std::shared_ptr<reference_branch> &parent,
                   DeltaType delta_from_parent)
      : iteration(parent->iteration + parent->size()), parent(parent),
        delta_to_parent(-delta_from_parent) {}

  // ref is the sequence of numbers at the center of this orbit
  void compute_terms(perturbation_orbit<OrbitType, DeltaType,
                                        HighPrecisionReferenceOrbit> &ref,
                     int max_iterations,
                     typename TermType::value_type max_delta_norm) {
    std::array<TermType, Terms> terms;
    terms[0] = 1.0;
    do {
      series.push_back({*ref, terms, ref});
      terms =
          HighPrecisionReferenceOrbit::calculation::delta_terms(*ref, terms);
      for (auto &t : terms) // TODO: Every 10 iterations
        t = normalize(t);
      ++ref;
    } while (iteration + size() < max_iterations && !escaped(*ref) &&
             max_delta_norm < series.back().maximum_delta_norm());
    final_entry = {*ref, terms, ref};
  }

  // The number of iterations covered by this branch
  int size() { return series.size(); }

  // On the first iteration, we have zero deviation from the epsilon
  // The z values in this series are the z values from THIS orbit
  // The terms in this series are the epsilon values from THIS orbit,
  // evaluated using deltas from THIS orbit.
  // This may be empty.
  std::vector<series_entry> series;

  // Our very last entry, which corresponds to the first iteration of any child
  // branches
  series_entry final_entry;

  // Delta is to the center of this branch
  int find_escape_iteration(DeltaType delta) const {
    // Look at the first term

    if (series.empty()) {
      if (parent)
        return parent->find_escape_iteration(delta + delta_to_parent);
      else
        return 0;
    }

    // Has the last term escaped
    if (!escaped(series.back().evaluate_z(delta)))
      return -1;

    // Look at the first item
    if (escaped(series.front().evaluate_z(delta)))
      return parent->find_escape_iteration(delta + delta_to_parent);

    int min = 0;
    int max = series.size() - 1;
    OrbitType z;

    while (max > min + 1) {
      auto mid = (min + max) / 2;
      assert(min >= 0 && mid < series.size());
      z = series[mid].evaluate_z(delta);
      if (escaped(z)) {
        max = mid - 1;
      } else {
        min = mid;
      }
    }
    return min + iteration;
  }
};

/*
This is an attempt to calculate the Mandelbrot set more efficiently. It
works by creating a recursive "tree" of reference orbits, and different
reference orbits can share the early terms in the reference orbit until it
is determined that the orbits need to branch.

We don't need to construct each tree explicitly, but rather follow the
outline of the tree using recursion.

Starting with a reference orbit, and a distance `epsilon` from that
reference orbit, we'll calculate a Taylor series for the relative epsilon
around the point `z+e`.

The "trick" is that we can use the divergence of the Taylor series to tell
us when to branch. We want to branch as late as possible as this allows us
to share more of the trees.

However, we don't need to store the terms of the Taylor series, we only need
to monitor if they diverge.

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

    central_delta: The position of this orbit relative to the central orbit.

    diagonal_size: The dimensions (width,height) of this region as a
        complex number.

central_orbit: A stored orbit generated from high-precision
        numbers. starts at {0,0} stop: A cancellation token

    iteration is initially 0.
    j is initially 0.
*/
template <Complex LowPrecisionType, Complex DeltaType, Complex TermType,
          unsigned long Terms, RandomAccessOrbit HighPrecisionReferenceOrbit,
          typename OutputFn>
void magic(int max_iterations, int x0, int y0, int w, int h,
           DeltaType diagonal_size, DeltaType central_delta,
           DeltaType delta_from_parent,
           const perturbation_orbit<LowPrecisionType, DeltaType,
                                    HighPrecisionReferenceOrbit> &parent_orbit,
           std::atomic<bool> &stop, OutputFn fn,
           const std::shared_ptr<
               reference_branch<LowPrecisionType, DeltaType, TermType, Terms,
                                HighPrecisionReferenceOrbit>> &parent) {

  using calculation = typename HighPrecisionReferenceOrbit::calculation;

  if (stop)
    return;
  const int Limit = 16;
  if (w <= Limit || h <= Limit) {

    // Let's dump the path
    for (auto p = parent; p; p = p->parent) {
      std::cout << (p->iteration + p->size()) << " "; //  << p->size();
      //
    }
    std::cout << std::endl;

    // The region we are rendering is very small, so we'll start calculating
    // each pixel, and be able to use the reference orbits we calculated along
    // the way

    // Look up the escape iteration instead using the branches
    // it = previous->find_escape_iteration(delta_from_previous);

    int iteration = parent ? parent->iteration + parent->size() : 0;

    auto pixel_delta =
        DeltaType{diagonal_size.real() / w, diagonal_size.imag() / h};
    auto half_diag =
        DeltaType{diagonal_size.real() * 0.5, diagonal_size.imag() * 0.5};
    auto top_left = central_delta - half_diag;

    for (int j = 0; j < h; ++j)
      for (int i = 0; i < w; ++i) {
        // auto delta_to_point = x;

        auto point_delta_from_parent =
            delta_from_parent - half_diag +
            DeltaType{pixel_delta.real() * i, pixel_delta.imag() * j};

        int it; // auto it =
                // parent->find_escape_iteration(point_delta_from_parent);
        if (true) {
          auto epsilon = parent->final_entry.epsilon(point_delta_from_parent);

          // To carry on the current orbit:
          auto orbit =
              parent_orbit.split_relative(point_delta_from_parent, epsilon);

          // To evaluate manually (slow):
          // perturbation_orbit<LowPrecisionType, DeltaType,
          //                   HighPrecisionReferenceOrbit>
          //    orbit(parent_orbit.reference,
          //          top_left + DeltaType{pixel_delta.real() * i,
          //                               pixel_delta.imag() * j});

          it = 0;

          for (; it < max_iterations && !escaped(*orbit); ++it)
            ++orbit;

          // auto it = 100 * fractals::norm(central_delta);

          // Iterate a small set of points
        }
        if (it < 0 || it >= max_iterations)
          it = 0;

        fn(x0 + i, y0 + j, it, iteration);
      }

    return;
  }

  // Create a new branch
  auto branch =
      parent ? std::make_shared<
                   reference_branch<LowPrecisionType, DeltaType, TermType,
                                    Terms, HighPrecisionReferenceOrbit>>(
                   parent, delta_from_parent)
             : std::make_shared<
                   reference_branch<LowPrecisionType, DeltaType, TermType,
                                    Terms, HighPrecisionReferenceOrbit>>();

  // Compute the terms of the new branch
  DeltaType epsilon_from_parent = {};

  if (parent) {
    epsilon_from_parent = parent->final_entry.epsilon(delta_from_parent);
  }

  // Pay attention: This bit is magic.
  // If we have one relative orbit, we can switch it to a different location
  // without recomputing all the terms from iteration 0.
  // We do need to know the new epsilon, which we computed using the Taylor
  // series.
  auto orbit =
      parent_orbit.split_relative(delta_from_parent, epsilon_from_parent);

  branch->compute_terms(orbit, max_iterations, fractals::norm(diagonal_size));

  // Branch step - split the current orbit into 4 parts
  auto half_size =
      DeltaType{diagonal_size.real() * 0.5, diagonal_size.imag() * 0.5};
  auto quarter_size = DeltaType{half_size.real() * 0.5, half_size.imag() * 0.5};

  auto lower_right_delta = quarter_size;
  auto upper_left_delta = -quarter_size;
  auto upper_right_delta =
      DeltaType(-upper_left_delta.real(), upper_left_delta.imag());
  auto lower_left_delta = -upper_right_delta;

  magic<LowPrecisionType, DeltaType, TermType, Terms>(
      max_iterations, x0, y0, w / 2, h / 2, half_size,
      central_delta + upper_left_delta, upper_left_delta, orbit, stop, fn,
      branch);
  magic<LowPrecisionType, DeltaType, TermType, Terms>(
      max_iterations, x0 + w / 2, y0, w - w / 2, h / 2, half_size,
      central_delta + upper_right_delta, upper_right_delta, orbit, stop, fn,
      branch);
  magic<LowPrecisionType, DeltaType, TermType, Terms>(
      max_iterations, x0, y0 + h / 2, w / 2, h - h / 2, half_size,
      central_delta + lower_left_delta, lower_left_delta, orbit, stop, fn,
      branch);
  magic<LowPrecisionType, DeltaType, TermType, Terms>(
      max_iterations, x0 + w / 2, y0 + h / 2, w - w / 2, h - h / 2, half_size,
      central_delta + lower_right_delta, lower_right_delta, orbit, stop, fn,
      branch);
}

template <Complex LowPrecisionType, Complex DeltaType, Complex TermType,
          unsigned long Terms, RandomAccessOrbit HighPrecisionReferenceOrbit,
          typename OutputFn>
void magic(int max_iterations, int x0, int y0, int w, int h,
           DeltaType diagonal_size,
           const HighPrecisionReferenceOrbit &central_orbit,
           std::atomic<bool> &stop, OutputFn fn) {

  perturbation_orbit<LowPrecisionType, DeltaType, HighPrecisionReferenceOrbit>
      orbit(central_orbit, {});

  magic<LowPrecisionType, DeltaType, TermType, Terms,
        HighPrecisionReferenceOrbit>(
      max_iterations, x0, y0, w, h, diagonal_size, {}, {}, orbit, stop, fn, {});
}
} // namespace mandelbrot