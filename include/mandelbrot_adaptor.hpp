#pragma once
#include "numbers/complex_number.hpp"

namespace mandelbrot {
// An "adaptor pattern" to transform the mandelbrot patterns
struct mandelbrot_adaptor {
  template <numbers::complex DeltaType, numbers::complex C>
  static const DeltaType &map_delta(const DeltaType &d, const C &c) {
    return d;
  }

  template <numbers::complex C> static const C &map(const C &p) { return p; }
};

struct mandeldrop_adaptor {
  template <numbers::complex C> static C map(const C &z) {
    auto n = numbers::inverse(numbers::norm(z));
    return {-imag(z) * n, -real(z) * n};
  }

  template <numbers::complex DeltaType, numbers::complex C>
  static DeltaType map_delta(const DeltaType &d, const C &c0) {
    // Note that we can calculate deltas using low precision (`double`) complex
    // numbers. We just need to rearrange our calculation of the delta to avoid
    // loss of precision.
    return DeltaType{0, 1} * d / numbers::number_cast<DeltaType>(c0 * (c0 + numbers::number_cast<C>(d)));
  }
};
} // namespace mandelbrot
