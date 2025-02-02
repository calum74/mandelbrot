#pragma once

#include "orbit.hpp"

namespace mandelbrot {

template <Complex DeltaType> DeltaType topleft(DeltaType radius) {
  return radius * -0.5;
}
template <Complex DeltaType> DeltaType topright(DeltaType radius) {
  return std::conj(radius) * 0.5;
}

template <Complex DeltaType> DeltaType bottomleft(DeltaType radius) {
  return -topright(radius);
}

template <Complex DeltaType> DeltaType bottomright(DeltaType radius) {
  return radius * 0.5;
}

// A type of perturbation orbit that supports
// "jump forward N steps"
template <Complex LowPrecisionType, Complex DeltaType, Complex TermType,
          RandomAccessOrbit ReferenceOrbit>
class orbit_branch {
public:
  /*
  TODO:
  Remove TermType and just use DeltaType
  Don't store final and base epsilon as they are in entries
  */

  using debug_orbit =
      perturbation_orbit<LowPrecisionType, DeltaType, ReferenceOrbit>;

  using reference_orbit_type = ReferenceOrbit;
  const ReferenceOrbit &reference_orbit;
  // Relative to the reference orbit
  DeltaType delta_from_reference, delta_from_parent;

  std::shared_ptr<const orbit_branch> parent;
  int base_iteration;

  struct entry {
    // The epsilon of this orbit against the reference orbit
    LowPrecisionType z;               // Of this orbit
    DeltaType epsilon_from_reference; // ?? Unclear if we need this
    DeltaType A, B;                   // ?? Not TermType
    int j;                            // ?? Maybe delete
    LowPrecisionType reference_z;
  };

  // ?? Do we even need to store all entries ?? For regions
  // that definitely don't escape, we can just skip this right?
  // We can probably throw this away and recompute it if needed.
  std::vector<entry> entries;

  using primary_orbit =
      perturbation_orbit<LowPrecisionType, DeltaType, ReferenceOrbit>;

  // A second order reference orbit
  // Problem: how do we rewind a perturbation orbit that's relative to another
  // one? using secondary_orbit =
  //    perturbation_orbit<LowPrecisionType, DeltaType, primary_orbit>;

  /*
    Gets the epsilon of point d relative to the reference orbit.

    Input:
    d - delta relative to this branch
    i - the absolute iteration number

    return - epsilon (dz) from the central orbit
  */
  DeltaType get_epsilon_1(int i, DeltaType delta_to_reference) const {

    if (i < base_iteration) {
      return parent->get_epsilon_1(i, delta_to_reference);
    } else if (parent) {
      return entries[i - base_iteration].epsilon_from_reference +
             entries[i - base_iteration].B *
                 (delta_to_reference - delta_from_reference);
    } else {

      // Optimization of previous case
      assert(base_iteration == 0);
      return entries[i].B * delta_to_reference;
    }
  }

  /*
    Same idea as get_epsilon_1, except that we'll recompute epsilon at each
    step.
  */
  DeltaType get_epsilon_2(int i, DeltaType delta_to_reference) const {

    if (i < base_iteration) {
      return parent->get_epsilon_2(i, delta_to_reference);
    } else if (parent) {
      auto e = parent->get_epsilon_2(base_iteration, delta_to_reference);
      return entries[i - base_iteration].A * e +
             entries[i - base_iteration].B *
                 (delta_to_reference - delta_from_reference);
    } else {
      assert(base_iteration == 0);
      return entries[i].B * delta_to_reference;
    }
  }

  // Fully precise version of the previous, for reference
  DeltaType get_epsilon_3(int i, DeltaType delta_to_reference) const {
    primary_orbit orbit(reference_orbit, delta_from_reference);
    while (i > 0) {
      i--;
      ++orbit;
    }
    return orbit.epsilon;
  }

  DeltaType get_epsilon(int i, DeltaType delta_to_reference) const {
    // Select which algorithm to use
    return get_epsilon_3(i, delta_to_reference);
  }

  typename DeltaType::value_type max_term_value(DeltaType radius) const {
    // return 100 * std::numeric_limits<typename
    // DeltaType::value_type>::epsilon();

    return std::numeric_limits<typename DeltaType::value_type>::epsilon() *
           std::numeric_limits<typename DeltaType::value_type>::epsilon() /
           fractals::norm(radius);
  }

  int get_escape_iterations(DeltaType d, int max_iterations,
                            int &skipped) const {

    auto i = size() + base_iteration;
    skipped = i;
    auto e = get_epsilon(i, d);

    // int j = i;
    // if (j >= reference_orbit.size() - 1)
    //   return 0;

    // Unfortunately it looks like we can't just reuse `j` from anywhere
    // Because the reference orbit should be against *our* z and we are doing a
    // double dereference at all times.

    int j = i;
    //    int j = entries.back().j; // [size()].j;
    auto z = e + entries.back().reference_z;

    if (!escaped(z)) { // !! Get an accurate z
      // Amazingly, we can restart the reference orbit at iteration 0 (j=0)
      // it doesn't actually matter what j is
      // The reference orbit is just there to avoid loss of precision
      perturbation_orbit<LowPrecisionType, DeltaType, ReferenceOrbit> orbit(
          reference_orbit, d, i, j, e); // Could be 0, z

      while (!escaped(*orbit) && i < max_iterations) {
        i++;
        ++orbit;
      }
      if (i == max_iterations)
        i = 0;
      return i;
    } else {
      // We'll need
    }

    return 0;

    //    return get_escape_iterations_naive(d, max_iterations);
  }

  // Keep this here for testing/debugging
  int get_escape_iterations_naive(DeltaType d, int max_iterations) const {

    perturbation_orbit<LowPrecisionType, DeltaType, ReferenceOrbit> orbit(
        reference_orbit, d);
    int i = 0;

    while (!escaped(*orbit) && i < max_iterations) {
      ++i;
      ++orbit;
    }
    if (i == max_iterations)
      i = 0;

    return i;
  }

  int size() const { return entries.size() - 1; }
  int j = 0; // ?? Private/delete?

  // Create the root, and populate it as far as it will go
  orbit_branch(const ReferenceOrbit &reference_orbit, DeltaType radius,
               int max_iterations, std::atomic<bool> &stop)
      : reference_orbit(reference_orbit), base_iteration(0) {
    delta_from_reference = 0;

    DeltaType epsilon = 0;
    auto max_B = max_term_value(radius);
    j = 0;

    DeltaType A{1}, B{0};

    LowPrecisionType z;

    do {
      z = reference_orbit[j];
      entries.push_back({z, epsilon, A, B, j, z});

      A = 2 * DeltaType{z} * A;
      B = 2 * DeltaType{z} * B + DeltaType{1, 0};
      j++;

      if (j >= reference_orbit.size() - 1 || escaped(reference_orbit[j]) ||
          fractals::norm(z) < fractals::norm(epsilon)) {
        epsilon = z;
        j = 0;
      }

    } while (!stop && !escaped(z) && j < reference_orbit.size() &&
             fractals::norm(B) <= max_B);
  }

  // Creates a branch, and populates it as far as it will go
  orbit_branch(const std::shared_ptr<const orbit_branch> &parent,
               DeltaType delta_from_parent, int max_iterations,
               std::atomic<bool> &stop)
      : reference_orbit(parent->reference_orbit), parent(parent),
        delta_from_parent(delta_from_parent),
        base_iteration(parent->base_iteration + parent->size()) {
    delta_from_reference = parent->delta_from_reference + delta_from_parent;
    // What is our epsilon at the b
    // base_epsilon =
    DeltaType epsilon =
        parent->get_epsilon(base_iteration, delta_from_reference);
    // j = parent->entries[parent->size()].j;
    j = parent->j;

    auto max_B = max_term_value(delta_from_parent);

    LowPrecisionType z;
    DeltaType A{1}, B{0};

    do {
      z = reference_orbit[j] + epsilon;
      entries.push_back({z, epsilon, A, B, j, reference_orbit[j]});

      A = 2 * DeltaType{z} * A;
      B = 2 * DeltaType{z} * B + DeltaType{1, 0};

      // Shouldn't need the final_epsilon squared term here??
      epsilon = 2 * z * epsilon + /* final_epsilon * final_epsilon + */
                delta_from_reference;
      j++;

      if (j == reference_orbit.size() - 1 || escaped(reference_orbit[j]) ||
          fractals::norm(z) < fractals::norm(epsilon)) {
        epsilon = z;
        j = 0;
      }

    } while (!stop && !escaped(z) && fractals::norm(B) <= max_B);
  }
};

template <typename Branch, Complex DeltaType, typename Fn>
void compute_tree(int x0, int y0, int x1, int y1,
                  const std::shared_ptr<Branch> &branch,
                  DeltaType branch_radius, int max_iterations,
                  std::atomic<bool> &stop, Fn fn) {
  if (stop)
    return;
  const int min = 32; // TODO: 8 probably

  if (x1 - x0 < min || y1 - y0 < min) {
    for (int j = y0; j < y1; ++j)
      for (int i = x0; i < x1; ++i) {
        DeltaType delta{(2.0 * double(i - x0) / double(x1 - x0) - 1.0) *
                            branch_radius.real(),
                        (2.0 * double(j - y0) / double(y1 - y0) - 1.0) *
                            branch_radius.imag()};
        int skipped;
        auto iterations = branch->get_escape_iterations(
            delta + branch->delta_from_reference, max_iterations, skipped);
        fn(i, j, iterations, skipped);
      }

    // Compute each pixel individually
  } else {
    // Call recursively with 4 more branches
    int mx = (x0 + x1) / 2;
    int my = (y0 + y1) / 2;
    auto new_radius = branch_radius * 0.5;

    compute_tree(x0, y0, mx, my,
                 std::make_shared<Branch>(branch, topleft(branch_radius),
                                          max_iterations, stop),
                 new_radius, max_iterations, stop, fn);

    compute_tree(mx, y0, x1, my,
                 std::make_shared<Branch>(branch, topright(branch_radius),
                                          max_iterations, stop),
                 new_radius, max_iterations, stop, fn);

    compute_tree(x0, my, mx, y1,
                 std::make_shared<Branch>(branch, bottomleft(branch_radius),
                                          max_iterations, stop),
                 new_radius, max_iterations, stop, fn);

    compute_tree(mx, my, x1, y1,
                 std::make_shared<Branch>(branch, bottomright(branch_radius),
                                          max_iterations, stop),
                 new_radius, max_iterations, stop, fn);
  }
}
} // namespace mandelbrot
