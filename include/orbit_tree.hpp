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

template <Complex DeltaType, Complex TermType>
std::array<TermType, 4> translate_terms(DeltaType delta,
                                        const std::array<TermType, 4> &terms) {
  // See orbit_tree.md for explanation
  auto delta_squared = fractals::square(delta);
  auto delta_cubed = delta_squared * delta;
  return {terms[0] + 2.0 * terms[1] * delta + 3.0 * terms[2] * delta_squared +
              4.0 * terms[3] * delta_cubed,
          terms[1] + 3.0 * terms[2] * delta + 6.0 * terms[3] * delta_squared,
          terms[2] + 4.0 * terms[3] * delta, terms[3]};
}

template <Complex LowPrecisionComplex, Complex DeltaType, Complex TermType,
          RandomAccessOrbit ReferenceOrbit, int Terms, int P1, int P2>
class orbit_branch {
public:
  using orbit_type =
      perturbation_orbit<LowPrecisionComplex, DeltaType, ReferenceOrbit>;
  using calculation = typename ReferenceOrbit::calculation;

  orbit_type orbit;

  int size() const { return entries.size() - 1; }

  // Create the root, and populate it as far as it will go
  orbit_branch(const ReferenceOrbit &reference_orbit, DeltaType radius,
               int max_iterations, std::atomic<bool> &stop) {

    base_iteration = 0;
    delta_from_reference_orbit = {0, 0};
    entries.push_back(entry()); // First iteration is zero
    orbit = orbit_type(reference_orbit, {});
    calculate_terms(radius, max_iterations, stop);
  }

  // Creates a branch, and populates it as far as it will go
  orbit_branch(const std::shared_ptr<const orbit_branch> &parent,
               DeltaType delta_from_parent, int max_iterations,
               std::atomic<bool> &stop)
      : parent(parent), delta_from_parent(delta_from_parent) {
    delta_from_reference_orbit =
        parent->delta_from_reference_orbit + delta_from_parent;
    base_iteration = parent->base_iteration + parent->size();
    auto radius_norm = fractals::norm(delta_from_parent);

    auto epsilon =
        evaluate_epsilon(delta_from_parent, parent->entries.back().terms);
    orbit = parent->orbit.split_relative(delta_from_parent, epsilon);
    entries.push_back({*orbit, translate_terms(delta_from_parent,
                                               parent->entries.back().terms)});
    calculate_terms(delta_from_parent, max_iterations, stop);
  }

  // TODO: Seed with previous value
  int get_escape_iterations(DeltaType delta, int max_iterations) const {

    // Case 1: The result is higher than the current branch
    if (!escaped(get_z(delta, size()))) {
      // The end hasn't escaped, so we need to keep iterating beyond the final
      // orbit Return the orbit with a delta
      auto orbit2 = orbit.split_relative(
          delta, evaluate_epsilon(delta, entries[size()].terms));
      while (orbit2.iteration() <= max_iterations && !escaped(*orbit2)) {
        ++orbit2;
      }
      return orbit2.iteration();
    }

    // Case 2: The result is lower than the current branch
    if (escaped(get_z(delta, 0))) {
      // Look in the parent
      if (parent) {
        return parent->get_escape_iterations(delta - delta_from_parent,
                                             max_iterations);
      }
      return 0;
    }

    // Case 3: The result is within the current branch
    int min = 0;
    int max = entries.size() - 1;
    while (min < max) {
      int mid = (min + max) / 2;
      if (escaped(get_z(delta, mid)))
        max = mid;
      else
        mid = mid;
    }
    return min;
  }

private:
  LowPrecisionComplex get_z(DeltaType delta, int i) const {
    return entries[i].z + evaluate_epsilon(delta, entries[i].terms);
  }

  void calculate_terms(DeltaType radius, int max_iterations,
                       std::atomic<bool> &stop) {
    auto radius_norm = fractals::convert<typename TermType::value_type>(
        fractals::norm(radius));
    do {
      ++orbit;
      entries.push_back(
          {*orbit, calculation::delta_terms(*orbit, entries.back().terms)});
      for (auto &t : entries.back().terms)
        t = normalize(t);
    } while (!stop && base_iteration + entries.size() <= max_iterations &&
             !escaped(*orbit) &&
             radius_norm < maximum_delta_norm<P1, P2>(entries.back().terms));
  }

  int base_iteration;
  std::shared_ptr<const orbit_branch> parent;
  DeltaType delta_from_parent, delta_from_reference_orbit;

  struct entry {
    LowPrecisionComplex z;
    std::array<TermType, Terms> terms;
  };

  std::vector<entry> entries;
};

template <typename Branch, Complex DeltaType, typename Fn>
void compute_tree(int x0, int y0, int x1, int y1,
                  const std::shared_ptr<Branch> &branch,
                  DeltaType branch_radius, int max_iterations,
                  std::atomic<bool> &stop, Fn fn) {
  const int min = 32;

  if (x1 - x0 < min || y1 - y0 < min) {
    // Compute each pixel individually
  } else {
    // Call recursively with 4 more branches
    int mx = (x0 + x1) / 2;
    int my = (y0 + y1) / 2;

    auto branch1 = std::make_shared<Branch>(branch, topleft(branch_radius),
                                            max_iterations, stop);
    compute_tree(x0, y0, mx, my, branch1, branch_radius * 0.5, max_iterations,
                 stop, fn);
  }
}
} // namespace mandelbrot
