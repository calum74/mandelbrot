#pragma once
#include "complex_number.hpp"

namespace mandelbrot {
// An "adaptor pattern" to transform the mandelbrot patterns
struct mandelbrot_adaptor {
  template <fractals::Complex DeltaType, fractals::Complex C>
  static const DeltaType &map_delta(const DeltaType &d, const C &c) {
    return d;
  }

  template <fractals::Complex C> static const C &map(const C &p) { return p; }
};

struct mandeldrop_adaptor {
  template <fractals::Complex C> static C map(const C &z) {
    auto n = fractals::inverse(fractals::norm(z));
    return {-imag(z) * n, -real(z) * n};
  }

  template <fractals::Complex DeltaType, fractals::Complex C>
  static DeltaType map_delta(const DeltaType &d, const C &c0) {
    // Note that we can calculate deltas using low precision (`double`) complex
    // numbers. We just need to rearrange our calculation of the delta to avoid
    // loss of precision.
    return DeltaType{0, 1} * d / fractals::convert<DeltaType>(c0 * (c0 + fractals::convert<C>(d)));
  }
};
} // namespace mandelbrot
