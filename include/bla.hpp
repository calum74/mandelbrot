#include "complex_number.hpp"
#include "orbit.hpp"

namespace mandelbrot {
  
template <typename SizeType>
SizeType max_norm_delta(SizeType normZ, SizeType normB)
{
  constexpr auto epsilon = std::numeric_limits<SizeType>::epsilon();
  constexpr auto e_squared = epsilon*epsilon;

  auto s1 = e_squared * normZ/normB;
  auto s2 = e_squared / (4 *  normB*normB);
  return std::min(s1, s2);
}

template <typename SizeType>
bool valid_approximation(SizeType normDelta, SizeType normZ, SizeType normB) {
  return normDelta < max_norm_delta(normZ, normB);
}

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
  int get(DeltaType delta, int max_iterations) {
    DeltaType epsilon = 0;
    int jZ = 0;
    DeltaType B = 0, B2=0;
    // 1) Find the top position
    auto norm_delta = norm(delta);

    // For now, just do a binary search
    int min = 0, max = stack.size() - 1;
    while (min + 5 < max) {
      int mid = (min + max) / 2;

      // Test the term at position mid to see if it still
      // has precision
      auto z = (*reference_orbit)[stack[mid].jZ];
      if(std::norm(stack[mid].B2) * norm(stack[mid].dc - delta) < 1e-7 * std::norm(stack[mid].B) )
//      if(norm(stack[mid].dc - delta) < stack[mid].maxNormDelta)
//      if (valid_approximation(norm(stack[mid].dc - delta), fractals::norm(z),
//                              norm(stack[mid].B)))
        min = mid;
      else
        max = mid;
    }

    last_jump = min;
    stack.resize(min);

    if (!stack.empty()) {
      B = stack.back().B;
      jZ = stack.back().jZ;
      // Use the delta_squared term??
      epsilon = B * (delta - stack.back().dc) + stack.back().dz;
    }

    // 2) Keep iterating until we escape

    using LowPrecisionType = typename Reference::value_type;
    LowPrecisionType z = (*reference_orbit)[jZ];

    typename DeltaType::value_type maxNormDelta;

    if(stack.empty())
      maxNormDelta = std::numeric_limits<typename DeltaType::value_type>::infinity();
    else
      maxNormDelta = stack.back().maxNormDelta;

    z = (*reference_orbit)[jZ] + LowPrecisionType(epsilon);
    do {
      epsilon =
          2 * (*reference_orbit)[jZ] * epsilon + epsilon * epsilon + delta;
      B2 = 2 * z * B2 + B*B;
      B = 2 * z * B + DeltaType{1, 0};

      maxNormDelta = std::min(maxNormDelta, max_norm_delta(fractals::norm(z), fractals::norm(B)));

      ++jZ;
      stack.push_back({B, B2, delta, epsilon, jZ, maxNormDelta});
      z = (*reference_orbit)[jZ] + LowPrecisionType(epsilon);

      // Zhuoran's device
      if (jZ >= reference_orbit->size() - 1 ||
          escaped((*reference_orbit)[jZ]) ||
          fractals::norm(z) < fractals::norm(epsilon)) {
        epsilon = z;
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
    typename DeltaType::value_type maxNormDelta;  // The size of validity
  };

  std::vector<entry> stack;
  int last_jump = 0;
};
} // namespace mandelbrot