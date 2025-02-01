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
  DeltaType base_epsilon, final_epsilon;
  DeltaType delta_from_reference, delta_from_parent;

  std::shared_ptr<const orbit_branch> parent;
  int base_iteration;

  struct entry {
    // The epsilon of this orbit against the reference orbit
    DeltaType epsilon_from_reference; // ?? Unclear if we need this
    DeltaType A, B;                   // ?? Not TermType
  };

  std::vector<entry> entries;

  // i is the absolute iteration number
  // d is the delta relative to this orbit
  // the return value is epsilon relative to the reference orbit
  DeltaType get_epsilon(int i, DeltaType d) const {
    assert(i >= 0 && i < base_iteration + entries.size());
    if (parent) {
      // Walk the tree from the root, computing the epsilon in the path
      // Get the original delta from the very root
      return entries[i - base_iteration].epsilon_from_reference +
             entries[i - base_iteration].A * get_base_epsilon(d) +
             entries[i - base_iteration].B * d;
    } else {
      assert(base_iteration == 0);
      return entries[i].B * d;
    }
  }

  // d must be in the "radius" of this orbit
  LowPrecisionType get_z(int i, DeltaType d) const {
    if (i < base_iteration) {
      return parent->get_z(i, d + delta_from_parent);
    }

    return get_epsilon(i, d) + reference_orbit[i]; // !! Should be entry.j
  }

  // If we had an orbit starting at d (relative to this orbit),
  // what epsilon should it have at the base of this branch?
  DeltaType get_base_epsilon(DeltaType d) const {
    return parent ? parent->get_epsilon(base_iteration, d + delta_from_parent)
                  : 0;
  }

  int get_escape_iterations(DeltaType d, int max_iterations) const {
    auto e = get_z(size(), d);

    perturbation_orbit<LowPrecisionType, DeltaType, ReferenceOrbit> orbit(
        reference_orbit, d + delta_from_reference);
    int i = 0;

    while (!escaped(*orbit) && i < max_iterations) {
      ++i;
      ++orbit;
    }
    if (i == max_iterations)
      i = 0;

    return i;
  }

  typename DeltaType::value_type max_term_value(DeltaType radius) const {
    return 1e-7 *
           std::numeric_limits<typename DeltaType::value_type>::epsilon() /
           fractals::norm(radius);
  }

  int size() const { return entries.size() - 1; }
  int j = 0;

  // Create the root, and populate it as far as it will go
  orbit_branch(const ReferenceOrbit &reference_orbit, DeltaType radius,
               int max_iterations, std::atomic<bool> &stop)
      : reference_orbit(reference_orbit), base_iteration(0) {
    base_epsilon = 0;

    auto max_B = max_term_value(radius);
    j = 0;

    DeltaType A{1}, B{0};

    auto z = reference_orbit[j];

    do {
      entries.push_back({{}, A, B});
      A = 2 * DeltaType{z} * A;
      B = 2 * DeltaType{z} * B + DeltaType{1, 0};
      z = reference_orbit[j];
      j++;

      if (j == reference_orbit.size() - 1 || escaped(reference_orbit[j]) ||
          fractals::norm(z) < fractals::norm(final_epsilon)) {
        final_epsilon = z;
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
    base_epsilon = parent->get_epsilon(base_iteration, delta_from_parent);
    final_epsilon = base_epsilon;
    j = parent->j;

#if 0 // Debug code only
    debug_orbit debug1{reference_orbit, delta_from_reference};
    for (int i = 0; i <= base_iteration; i++)
      ++debug1;
    debug_orbit debug2{reference_orbit, delta_from_reference, base_iteration, j,
                       final_epsilon};

    // debug1 and debug2 should now be identical to this orbit
#endif
    auto max_B = max_term_value(delta_from_parent);

    LowPrecisionType z = reference_orbit[j] + final_epsilon;
    DeltaType A{1}, B{0};

    do {
      // assert(final_epsilon == debug1.epsilon);
      // assert(final_epsilon == debug2.epsilon);
      //++debug1;
      //++debug2;
      entries.push_back({final_epsilon, A, B});
      A = 2 * DeltaType{z} * A;
      B = 2 * DeltaType{z} * B + DeltaType{1, 0};

      // Shouldn't need the final_epsilon squared term here??
      final_epsilon =
          2 * z * final_epsilon + /* final_epsilon * final_epsilon + */
          delta_from_reference;
      j++;
      z = reference_orbit[j] + final_epsilon;

      if (j == reference_orbit.size() - 1 || escaped(reference_orbit[j]) ||
          fractals::norm(z) < fractals::norm(final_epsilon)) {
        final_epsilon = z;
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
        fn(i, j, branch->get_escape_iterations(delta, max_iterations));
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
