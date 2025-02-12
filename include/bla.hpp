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
A bilinear orbit with a fixed step size.
At each step, you can jump forward at most `step_size` iterations.


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
    TermType A{1}, A2{0}, B{0}, B2{0}, C{0};
    // 1) Find the top position

    // The entry at "n" allows us to skip forward to "n" from a previous step
    int n = step_size;
    int min = 0;

    // dz is always relative to the current orbit
    while (n < stack.size()) {
      TermType local_dc = dc - stack[n].dc;
      using term_real = typename TermType::value_type;
      auto & entry = stack[n];

      // TODO: Use termtype comparisons but unfortunately there's a bug in <=
      // double b1 = convert<double>(fractals::norm(stack[n].B2 * TermType(local_dc)));
      // double b2 = convert<double>(term_real(1e-7) * fractals::norm(stack[n].B));
      // double a1 = convert<double>(fractals::norm(stack[n].A2 * TermType(dz)));
      // double a2 = convert<double>(term_real(1e-7) * fractals::norm(stack[n].A));

      auto e1 = entry.A * TermType{dz} + entry.B * local_dc;
      auto e2 = entry.A2 * local_dc * local_dc + entry.B2 * TermType(dz * dz) + local_dc * TermType(dz) * entry.C;
      auto x = std::numeric_limits<double>::epsilon();

      if(convert<double>(fractals::norm(e2)) <= x * convert<double>(fractals::norm(e1))) {
        // fractals::norm(e2) < 1e-7 * fractals::norm(e1)) {
        // !! Check the validity ?? What about C
        // if (b1<=b2 && a1<=a2) {
        dz = convert<DeltaType>(entry.A * TermType{dz} + entry.B * local_dc);  // !! = e1
        min = n;
        n += step_size;

        // !! Be careful with orbit resetting I think
      } else {
        // std::cout << e1 << ' ' << e2 << std::endl;
        // std::cout << a1 << " " << a2 << " " << b1 << " " << b2 << std::endl;
        break;
      }
    }

    // TODO: We could try to skip forward further in the final segment

    last_jump = min;
    stack.resize(min);

    // Our very last stack entry allows us to continue the sequence
    if (!stack.empty()) {
      A = stack.back().A;
      A2 = stack.back().A2;
      B = stack.back().B;
      B2 = stack.back().B2;
      C = stack.back().C;
      jZ = stack.back().jZ;

      // Note that by adding stack.back().dz
      // we have translated the epsilon to be relative to the central orbit

      // Now, dz is relative to the high-precision orbit
      dz = convert<DeltaType>(A * TermType(dz) +
                              B * TermType(dc - stack.back().dc)) +
           stack.back().dz;
    }

    // 2) Keep iterating until we escape or reach the iteration limit

    using LowPrecisionType = typename Reference::value_type;
    LowPrecisionType z = (*reference_orbit)[jZ] + convert<LowPrecisionType>(dz);
    do {

      // z = (*reference_orbit)[jZ] + LowPrecisionType(dz);

      // ?? Where does this go
      stack.push_back({A, A2, B, B2, C, dc, dz, jZ});

      dz = 2 * (*reference_orbit)[jZ] * dz + dz * dz + dc;

      ++jZ;


      if (stack.size() % step_size == 0) {
        // Reset the bivariate coefficients
        // Including quadratic terms
        A = 1;
        A2 = B = B2 = C = 0;
      } 
      // else  // ?? Do we need the else?
       {
        // !! Note the order of these
        TermType zz = 2 * z; // ?? Do we just want (*reference_orbit)[jZ] here?
        C = zz * C + 2 * A * B;
        A2 = zz * A2 + A * A;
        B2 = zz * B2 + B * B;
        A = zz * A;
        B = zz * B + TermType{1, 0};
      }

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
    TermType A, A2, B, B2, C;
    DeltaType dc, dz;
    int jZ; // Zhouran's j. It's probably implicit from the entry position but
            // save it for now.
  };

  std::vector<entry> stack;
  int last_jump = 0;
};
} // namespace mandelbrot