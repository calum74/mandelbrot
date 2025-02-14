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
             fractals::norm(p.first);
}

template <typename DeltaType, typename TermType> class linear_step {
public:
  linear_step() = default;

  linear_step(DeltaType dc)
      : A{1}, A2{}, B{}, B2{}, C{}, dc{dc}, dz{}, jZ(0), debug_from(0), debug_to(0) {}

  linear_step restart() const {
    return {TermType{1}, TermType{0}, TermType{0}, TermType{0}, TermType{0},
            dc,          dz,          jZ,          debug_to,  debug_to};
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
  std::optional<DeltaType> jump_dz(DeltaType dz1, DeltaType dc1, int expected_from, int expected_to) const {
    assert(debug_from == expected_from);
    assert(debug_to == expected_to);
    auto l = get_local_dz(dz1, dc1 - dc);
    return is_valid(l)
               ? std::optional<DeltaType>{convert<DeltaType>(l.first)}
               : std::nullopt;
  }

  // Note that new_dc and new_dz are relative to the reference orbit
  linear_step move_to(DeltaType new_dz, DeltaType new_dc) const {
    /* I think only the linear terms are valid to move, we can also move the
    higher order terms as they are only used for error estimation. More
    investigation needed.
    */

    return {A, A2, B, B2, C, new_dc, new_dz, jZ, debug_from, debug_to};
  }

  template <RandomAccessOrbit Reference>
  auto get_z(const Reference &reference) const {
    return reference[jZ] + dz;
  }

private:
public:  // Temporary
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
At each step, you can jump forward at most `step_size` iterations.

`stack[i]` contains terms needed to jump forward from iteration `i%step_size` to
iteration `i`. If `i%step_size == 0`, then this jumps forward from
`i-step_size`.

Iteration 0 is unused since you'll never need to jump forward to 0.

To jump, calculate `dz = e.A*dz + e.B*(dc-e.dc)`, after first checking the
higher order terms (A2,B2,C) for divergence.
*/
template <Complex DeltaType, Complex TermType, RandomAccessOrbit Reference>
class bilinear_orbit {
public:
  bilinear_orbit() : reference_orbit() {}
  bilinear_orbit(const Reference &r) : reference_orbit(&r) {}

  static constexpr int step_size = 50;
  using entry = linear_step<DeltaType, TermType>;

  // The delta dc is to the reference orbit
  int get(const DeltaType dc, int max_iterations) {

    // !! Don't reuse dz for both orbits
    DeltaType dz = 0; // dz from the reference orbit
    int jZ = 0;
    // 1) Find the top position

    // The entry at "n" allows us to skip forward to "n" from a previous step
    int n = step_size;
    int min = 0;

    // dz is always relative to the current orbit
    while (n < stack.size()) {
      // there's some mismatch between the jump_dz and the computed dz
      // off by one I think!

      auto j = stack[n].jump_dz(dz, dc, n-step_size, n);

      if (j) {
        dz = *j;
        min = n;
        n += step_size;
        if (n == 3*step_size) break; // !! Temporary
      } else {
        break;
      }
    }

    // TODO: We could try to skip forward further in the final segment
    last_jump = min;

    int debug_from, debug_to;

    if (stack.empty()) {
      stack.push_back(entry(dc));
    } else {
      stack.resize(min + 1);
      dz += stack.back().dz;
    }

    // 2) Keep iterating until we escape or reach the iteration limit

    entry e = stack.back().move_to(dz, dc);

    do {
      if (stack.size() % step_size == 1) {
        e = e.restart();
      }
      e = e.next_iteration(*reference_orbit);
      stack.push_back(e);

    } while (!escaped(stack.back().get_z(*reference_orbit)) &&
             stack.size() < max_iterations);

    if (stack.size() == max_iterations)
      return 0;

    return stack.size();
  }

  int get_skipped_iterations() const { return last_jump; }

private:
  const Reference *reference_orbit;

  std::vector<entry> stack;
  int last_jump = 0;
};
} // namespace mandelbrot