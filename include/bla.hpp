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
  TermType A = 1, A2 = 0, B = 0, B2 = 0, C = 0;
};

template <typename TermType>
jump_terms<TermType> step_terms(const jump_terms<TermType> &terms, TermType z) {
  return {2 * z * terms.A, 2 * z * terms.A2 + terms.A * terms.A,
          2 * z * terms.B + 1, 2 * z * terms.B2 + terms.B * terms.B,
          2 * z * terms.C + 2 * terms.A * terms.B};
}

template <typename DeltaType, typename TermType> class jump_step {
public:
  int n; // The iteration we are jumping from
  int m; // The iteration we are jumping to
  jump_terms<TermType> terms;
  DeltaType d_ik, e_nik, e_mik;
};

template <typename DeltaType, typename TermType> class linear_step {
public:
  linear_step() = default;

  linear_step(DeltaType dc)
      : A{1}, A2{}, B{}, B2{}, C{}, dc{dc}, dz{}, jZ(0), debug_from(0),
        debug_to(0) {}

  linear_step restart() const {
    return {TermType{1}, TermType{0}, TermType{0}, TermType{0}, TermType{0},
            dc,          dz,          jZ,          debug_to,    debug_to};
  }

  template <RandomAccessOrbit Reference>
  linear_step next_iteration(const Reference &reference) const {
    auto z = reference[jZ];
    auto orbit_z = z + dz;
    auto dz_next = 2 * z * dz + dz * dz + dc;

    int jZ_next = jZ + 1;

    if (jZ_next >= reference.size() - 1 || escaped(z) ||
        fractals::norm(z) < fractals::norm(dz)) {
      dz_next =
          reference[jZ_next] + convert<typename Reference::value_type>(dz_next);
      jZ_next = 0;
    }

    TermType two_orbit_z = 2 * orbit_z;

    return {two_orbit_z * A,
            two_orbit_z * A2 + A * A,
            two_orbit_z * B + TermType{1, 0},
            two_orbit_z * B2 + B * B,
            two_orbit_z * C + 2 * A * B,
            dc,
            dz_next,
            jZ_next,
            debug_from,
            debug_to + 1};
  }

  // Returns a dz (relative to the reference orbit)
  // If the precision is invalid, returns nothing
  // ?? Return complex NaN?
  std::optional<DeltaType> jump_dz(DeltaType dz1, DeltaType dc1,
                                   int expected_from, int expected_to) const {
    assert(debug_from == expected_from);
    assert(debug_to == expected_to);
    auto l = get_local_dz(dz1, dc1 - dc);
    using T = double;

    if (convert<T>(fractals::norm(B2)) * fractals::norm(dc1 - dc) >
        1e-15 * convert<T>(fractals::norm(B)))
      return std::nullopt;

    if (convert<T>(fractals::norm(A2)) * fractals::norm(dz1) >
        1e-15 * convert<T>(fractals::norm(A)))
      return std::nullopt;

    // !! Test for escape
    return is_valid(l) ? std::optional<DeltaType>{convert<DeltaType>(l.first)}
                       : std::nullopt;
  }

  template <RandomAccessOrbit Reference>
  auto get_z(const Reference &reference) const {
    return reference[jZ] + dz;
  }

private:
public: // Temporary
  linear_step(TermType a, TermType a2, TermType b, TermType b2, TermType c,
              DeltaType dc, DeltaType dz, int jZ, int from, int to)
      : A(fractals::normalize(a)), A2(fractals::normalize(a2)),
        B(fractals::normalize(b)), B2(fractals::normalize(b2)),
        C(fractals::normalize(c)), dc(dc), dz(dz), jZ(jZ), debug_from(from),
        debug_to(to) {}

  TermType A, A2, B, B2, C;
  DeltaType dc, dz;
  int jZ;                   // Index into the reference orbit
  int debug_from, debug_to; // Not needed but useful for checking

  // Gets the new dz, relative to the orbit from our dc.
  // You would need to add dz to get the orbit relative to the reference orbit.
  // !! TermType (except that it doesn't work!)
  std::pair<TermType, TermType> get_local_dz(TermType dz1, TermType dc1) const {
    return {A * dz1 + B * dc1, A2 * dz1 * dz1 + B2 * dc1 * dc1 + C * dc1 * dz1};
  }
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

    int n=0;
    DeltaType e_njk = 0;
    int reference_index = 0;

    using LowPrecisionType = typename Reference::value_type;
    LowPrecisionType z_nj = (*reference_orbit)[reference_index] + e_njk;

    while(!escaped(z_nj) && n < max_iterations)
    {
      // Perform once step
      auto z_nk = (*reference_orbit)[reference_index];

      // Next iteration
      e_njk = 2 * z_nk * e_njk + e_njk*e_njk + d_jk;

      ++reference_index;
      ++n;

      z_nj = (*reference_orbit)[reference_index] + e_njk;

      if(reference_index >= reference_orbit->size()-1)
      {
        reference_index = 0;
        e_njk = z_nj;
      }
    }

    if(n == max_iterations) return 0;

    return n;
  }

  int get_skipped_iterations() const { return last_jump; }

private:
  const Reference *reference_orbit;

  std::vector<step> steps;

  int last_jump = 0;
};
} // namespace mandelbrot