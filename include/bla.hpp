#include "complex_number.hpp"
#include "orbit.hpp"

#include <iostream>

namespace mandelbrot {

/*
A linear approximation orbit.

Another experiment.

Simpler than bilinear, because this orbit only ever starts from iteration 0
so the "dz" term is always 0.

It may seem a bit of a cowboy move to reuse the lower terms of an unrelated
orbit, but the key insight is that for "close" orbits, the (bi)linear terms are
provably identical (within epsilon). We'll chop off the top of the orbit and
recalculate the terms that might have changed.

This orbit is more expensive than a regular perturbation orbit because we
construct a vector of the entire orbit, but the hope is that we can reuse most
of the orbit for adjacent orbits.
*/
template <Complex DeltaType, RandomAccessOrbit Reference> class linear_orbit {
public:
  linear_orbit() : reference_orbit() {}
  linear_orbit(const Reference &r) : reference_orbit(&r) {}

  // The delta is to the reference orbit
  int get(DeltaType dc, int max_iterations) {
    DeltaType dz = 0;
    int jZ = 0;
    DeltaType B = 0, B2 = 0;
    // 1) Find the top position
    auto norm_dc = norm(dc);

    // For now, just do a binary search
    int min = 0, max = stack.size() - 1;
    while (min + 5 < max) {
      int mid = (min + max) / 2;

      // Test the term at position mid to see if it still
      // has precision
      auto z = (*reference_orbit)[stack[mid].jZ];
      if (std::norm(stack[mid].B2) * norm(stack[mid].dc - dc) <
          1e-7 * std::norm(stack[mid].B))
        min = mid;
      else
        max = mid;
    }

    last_jump = min;
    stack.resize(min);

    if (!stack.empty()) {
      B = stack.back().B;
      B2 = stack.back().B2;
      jZ = stack.back().jZ;
      // Use the delta_squared term??
      dz = B * (dc - stack.back().dc) + stack.back().dz;
    }

    // 2) Keep iterating until we escape

    using LowPrecisionType = typename Reference::value_type;
    LowPrecisionType z = (*reference_orbit)[jZ];

    z = (*reference_orbit)[jZ] + LowPrecisionType(dz);
    do {
      dz = 2 * (*reference_orbit)[jZ] * dz + dz * dz + dc;
      B2 = 2 * z * B2 + B * B;
      B = 2 * z * B + DeltaType{1, 0};

      ++jZ;
      stack.push_back({B, B2, dc, dz, jZ});
      z = (*reference_orbit)[jZ] + LowPrecisionType(dz);

      // Zhuoran's device
      if (jZ >= reference_orbit->size() - 1 ||
          escaped((*reference_orbit)[jZ]) ||
          fractals::norm(z) < fractals::norm(dz)) {
        dz = z;
        jZ = 0;
      }

    } while (!escaped(z) && stack.size() < max_iterations);

    if (stack.size() == max_iterations)
      return 0;

    return stack.size();
  }

  int get_skipped_iterations() const { return last_jump; }

private:
  const Reference *reference_orbit;

  struct entry {
    DeltaType B, B2, dc, dz;
    int jZ; // Zhouran's j. It's probably implicit from the entry position but
            // save it for now.
  };

  std::vector<entry> stack;
  int last_jump = 0;
};

template <typename DeltaType> bool is_valid(std::pair<DeltaType, DeltaType> p) {
  return fractals::norm(p.second) <=
         std::numeric_limits<typename DeltaType::value_type>::epsilon() *
             std::numeric_limits<typename DeltaType::value_type>::epsilon() *
             fractals::norm(p.first);
}

template <typename TermType> struct jump_terms {
  TermType A = {1}, A2 = {0}, B = {0}, B2 = {0}, C = {0};
};

template <typename TermType>
jump_terms<TermType> step_terms(const jump_terms<TermType> &terms, TermType z) {
  return {2 * z * terms.A, 2 * z * terms.A2 + terms.A * terms.A,
          2 * z * terms.B + TermType{1}, 2 * z * terms.B2 + terms.B * terms.B,
          2 * z * terms.C + 2 * terms.A * terms.B};
}

template <typename DeltaType, typename TermType> class jump_step {
public:
  int n; // The iteration we are jumping from
  int m; // The iteration we are jumping to
  jump_terms<TermType> terms;
  DeltaType d_ik, e_nik, e_mik;
};



/*
A bivariate orbit with a fixed step size.
*/
template <Complex DeltaType, Complex TermType, RandomAccessOrbit Reference>
class bilinear_orbit {
public:
  bilinear_orbit() : reference_orbit() {}
  bilinear_orbit(const Reference &r) : reference_orbit(&r) {}

  static constexpr int step_size = 50;
  using step = jump_step<DeltaType, TermType>;

  // k is the reference orbit
  // j is the current orbit
  int get(const DeltaType d_jk, int max_iterations) {

    int n = 0;
    DeltaType e_njk = 0;
    int reference_index = 0;

    using LowPrecisionType = typename Reference::value_type;
    LowPrecisionType z_nj = (*reference_orbit)[reference_index] + e_njk;

    step new_step;

    while (!escaped(z_nj) && n < max_iterations) {
      if (n % step_size == 0) {
        if (n > 0) {
          int index = n / step_size;
          new_step.e_mik = e_njk;
          new_step.m = n;

          // Overwrite/push current term
          if (index < steps.size())
            steps[index] = new_step;
          else
            steps.push_back(new_step);
        }

        // See how far we can skip forward
        // TODO

        // Reset the terms for manual calculation
        new_step = {.n=n, .d_ik = d_jk, .e_nik = e_njk};
      }

      // Perform one step

      new_step.terms = step_terms(new_step.terms, TermType(z_nj));

      auto z_nk = (*reference_orbit)[reference_index];

      // Next iteration

      ++reference_index;
      ++n;

      e_njk = 2 * z_nk * e_njk + e_njk * e_njk + d_jk;
      z_nj = (*reference_orbit)[reference_index] + e_njk;

      if (reference_index >= reference_orbit->size() - 1) {
        reference_index = 0;
        e_njk = z_nj;
      }
    }

    if (n == max_iterations)
      return 0;

    return n;
  }

  int get_skipped_iterations() const { return skipped_steps; }

private:
  const Reference *reference_orbit;

  std::vector<step> steps;

  int skipped_steps = 0;
};
} // namespace mandelbrot