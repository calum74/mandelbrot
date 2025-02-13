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

  static constexpr int step_size = 20;

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
      TermType local_dc = dc - stack[n].dc;
      using term_real = typename TermType::value_type;
      auto &entry = stack[n];

      assert(entry.debug_from == min);
      assert(entry.debug_to == n);

      auto new_dz = entry.A * TermType{dz} + entry.B * local_dc;
      auto residual = entry.A2 * local_dc * local_dc +
                      entry.B2 * TermType(dz * dz) +
                      entry.C * local_dc * TermType(dz);
      auto E = std::numeric_limits<double>::epsilon();

      // !! convert needed since <= does not work properly for 0
      if (convert<double>(fractals::norm(residual)) <=
          E * convert<double>(fractals::norm(new_dz))) {
        dz = convert<DeltaType>(new_dz);
        min = n;
        n += step_size;

        break; // Do not check in

        // !! Be careful with orbit resetting I think
      } else {
        // std::cout << e1 << ' ' << e2 << std::endl;
        // std::cout << a1 << " " << a2 << " " << b1 << " " << b2 << std::endl;
        break;
      }
    }

    // TODO: We could try to skip forward further in the final segment
    last_jump = min;

    int debug_from, debug_to;
    min = 0; // Test

    if (stack.empty()) {
      debug_from = 0;
      debug_to = 0;
      TermType A{1}, A2{0}, B{0}, B2{0}, C{0};
      stack.push_back({A, A2, B, B2, C, dc, dz, jZ, 0, 0});
    } else {
      // Our very last stack entry allows us to continue the sequence
      stack.resize(min + 1);
      jZ = stack.back().jZ;
      debug_from = stack.back().debug_from;
      debug_to = stack.back().debug_to;

      // Note that by adding stack.back().dz
      // we have translated the epsilon to be relative to the central orbit

      // Now, dz is relative to the high-precision orbit
      dz = convert<DeltaType>(stack.back().A * TermType(dz) +
                              stack.back().B * TermType(dc - stack.back().dc)) +
           stack.back().dz;
    }

    // 2) Keep iterating until we escape or reach the iteration limit

    using LowPrecisionType = typename Reference::value_type;
    LowPrecisionType z = (*reference_orbit)[jZ] + convert<LowPrecisionType>(dz);
    do {

      // Our current iteration
      // Iteration 0 is at zero.
      int i = stack.size() - 1;

      // The distance of the current orbit to the reference orbit at the current
      // iteration

      auto dz_i = dz;
      auto &previous_entry = stack.back();

      TermType A_in, A2_in, B_in, B2_in, C_in;
      if (i % step_size == 0) {
        A_in = 1;
        A2_in = B_in = B2_in = C_in = 0;
        debug_from = i;
      } else {
        A_in = previous_entry.A;
        A2_in = previous_entry.A2;
        B_in = previous_entry.B;
        B2_in = previous_entry.B2;
        C_in = previous_entry.C;
      }

      auto z_i = (*reference_orbit)[jZ];
      auto dz_i1 = 2 * z_i * dz_i + dz_i * dz_i + dc;

      auto local_z_i = z_i + dz_i;

      TermType two_z_i = 2 * local_z_i;

      auto A_in1 = two_z_i * A_in;
      auto A2_in1 = two_z_i * A2_in + A_in * A_in;
      auto B_in1 = two_z_i * B_in + TermType{1, 0};
      auto B2_in1 = two_z_i * B2_in + B_in * B_in;
      auto C_in1 = two_z_i * C_in + 2 * A_in * B_in;

      stack.push_back({A_in1, A2_in1, B_in1, B2_in1, C_in1, dc, dz_i1, jZ,
                       debug_from, ++debug_to});

      // Next iteration
      ++jZ;

      dz = dz_i1;
      // Zhuoran's device
      z = (*reference_orbit)[jZ] + LowPrecisionType(dz);
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
    entry(TermType a, TermType a2, TermType b, TermType b2, TermType c,
          DeltaType dc, DeltaType dz, int jZ, int from, int to)
        : A(fractals::normalize(a)), A2(fractals::normalize(a2)),
          B(fractals::normalize(b)), B2(fractals::normalize(b2)),
          C(fractals::normalize(c)), dc(dc), dz(dz), jZ(jZ), debug_from(from),
          debug_to(to) {}
    entry() = default;
    entry(const entry &) = default;
    TermType A, A2, B, B2, C;
    DeltaType dc, dz;
    int jZ; // Zhouran's j. It's probably implicit from the entry position but
            // save it for now.
    int debug_from, debug_to; // Not needed but useful for checking
  };

  std::vector<entry> stack;
  int last_jump = 0;
};
} // namespace mandelbrot