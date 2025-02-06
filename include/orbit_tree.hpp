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
    DeltaType epsilon_from_reference; // dz from the primary reference orbit
    DeltaType A, B; // Bilinear coefficients from the secondary orbit
    int j; // The index of the reference orbit (wrapped using Zhouran's device)
  };

  // ?? Do we even need to store all entries ?? For regions
  // that definitely don't escape, we can just skip this right?
  // We can probably throw this away and recompute it if needed.
  std::vector<entry> entries;

  using primary_orbit =
      perturbation_orbit<LowPrecisionType, DeltaType, ReferenceOrbit>;

  /*
    Same idea as get_epsilon_1, except that we'll recompute epsilon at each
    step.
  */
  std::pair<DeltaType, int> get_epsilon_2(int i,
                                          DeltaType delta_reference) const {

    if (i < base_iteration) {
      return parent->get_epsilon_2(i, delta_reference);
    } else if (parent) {

      if (size() == 0)
        return parent->get_epsilon_2(i, delta_reference);

      auto &entry = entries[i - base_iteration];
      auto epsilon = parent->get_epsilon_2(base_iteration, delta_reference);
      return std::make_pair(
          entry.A * (epsilon.first - entries[0].epsilon_from_reference) +
              entry.B * (delta_reference - delta_from_reference) +
              entry.epsilon_from_reference,
          entry.j);
    } else {
      assert(base_iteration == 0);
      auto &entry = entries[i - base_iteration];
      return std::make_pair(
          entry.B * delta_reference + entry.epsilon_from_reference, entry.j);
    }
  }

  // Fully precise version of the previous, for reference
  DeltaType get_epsilon_3(int i, DeltaType delta_to_reference) const {
    primary_orbit orbit(reference_orbit, delta_to_reference);
    while (i > 0) {
      i--;
      ++orbit;
    }
    return std::make_pair(orbit.epsilon, i);
  }

  std::pair<DeltaType, int> get_epsilon(int i,
                                        DeltaType delta_reference) const {
    // Select which algorithm to use
    return get_epsilon_2(i, delta_reference);
  }

  int get_escape_iterations(DeltaType d, int max_iterations,
                            int &skipped) const {

    auto i = size() + base_iteration;
    skipped = i;
    auto e = get_epsilon(i, d);

    int j = e.second;
    auto z = e.first + reference_orbit[j];

    if (!escaped(z)) {
      perturbation_orbit<LowPrecisionType, DeltaType, ReferenceOrbit> orbit(
          reference_orbit, d, i, e.second, e.first);

      while (!escaped(*orbit) && i < max_iterations) {
        i++;
        ++orbit;
      }
      if (i == max_iterations)
        i = 0;
      return i;
    } else {
      // Conclusion: We probably don't need to store all the vectors after all
      return 0;
      int min = 0;
      int max = i;
      while (min < max) {
        int mid = (min + max) / 2;
        auto e = get_epsilon(mid, d);
        auto z = e.first + reference_orbit[e.second];
        if (escaped(z))
          max = mid;
        else
          min = mid + 1;
      }
      return min;
    }

    return 0;
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

  // Create the root, and populate it as far as it will go
  orbit_branch(const ReferenceOrbit &reference_orbit, DeltaType radius,
               int max_iterations, std::atomic<bool> &stop)
      : reference_orbit(reference_orbit), base_iteration(0) {
    delta_from_reference = 0;

    DeltaType epsilon = 0;
    int j = 0;

    extend_series(j, epsilon, fractals::norm(radius), max_iterations, stop);
  }

  using delta_size = typename DeltaType::value_type;

  void extend_series(int j, DeltaType epsilon, auto norm_delta,
                     int max_iterations, std::atomic<bool> &stop) {

    DeltaType A{1}, B{0};

    DeltaType local_epsilon{0}, local_A{1}, local_B{0};
    LowPrecisionType local_z;

    LowPrecisionType z;

    // Tweak to add extra inaccuracy
    constexpr auto system_epsilon = // 1e-6;
        std::numeric_limits<typename DeltaType::value_type>::epsilon();
    constexpr auto system_epsilon_squared = system_epsilon * system_epsilon;

    using size_type = typename DeltaType::value_type;
    delta_size max_norm_A, max_norm_B;
    delta_size max_parent_epsilon = parent ? parent->max_epsilon : size_type{0};
    // auto modulus_delta = std::sqrt(norm_delta);

    do {
      z = reference_orbit[j];
      local_z = z + epsilon;

      entries.push_back({epsilon, A, B, j});

      // The A and B terms are for linearization around this orbit
      // Meanwhile, we're also computing our own orbit using perturbation theory
      A = 2 * DeltaType{local_z} * A;
      B = 2 * DeltaType{local_z} * B + DeltaType{1, 0};

      epsilon = 2 * z * epsilon +
                epsilon * epsilon + /* ?? Delete this term surely final_epsilon
                                     * final_epsilon + */
                delta_from_reference;
      j++;

      auto norm_z = fractals::norm(local_z);

      if (norm_z == 0)
        norm_z = 1;

      // Zhuoran's device
      if (j >= reference_orbit.size() - 1 || escaped(reference_orbit[j]) ||
          norm_z < fractals::norm(epsilon)) {
        epsilon = z;
        j = 0;
      }

      max_norm_A = 2 * system_epsilon_squared * norm_z /
                   (max_parent_epsilon * max_parent_epsilon);

      max_norm_B = system_epsilon_squared * norm_z / norm_delta;

    } while (!stop && !escaped(local_z) &&
             entries.size() + base_iteration < max_iterations &&
             norm(A) < max_norm_A && norm(B) < max_norm_B);

    max_epsilon = std::sqrt(norm_delta) * std::sqrt(fractals::norm(B));

    if (!parent || !parent->parent) {
      std::cout << "||∂|| = " << norm_delta;
      std::cout << ", Size = " << size() << ", Norm A = " << norm(A)
                << ", max norm A = " << max_norm_A << ", Norm B = " << norm(B)
                << ", max norm B = " << max_norm_B << std::endl;
    }
    // std::cout << max_epsilon << std::endl;

    if (parent)
      max_epsilon += std::sqrt(norm(A)) * parent->max_epsilon;
  }

  // Note our calculation for this is bogus !!
  delta_size max_epsilon;

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
        parent->get_epsilon(base_iteration, delta_from_reference).first;
    // The first iteration of this branch is the size() iteration from the
    // parent branch
    int j = parent->entries[parent->size()].j;

    extend_series(j, epsilon, fractals::norm(delta_from_parent), max_iterations,
                  stop);
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
        int skipped = 0;
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
