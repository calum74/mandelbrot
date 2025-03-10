#pragma once
#include <cstdint>

namespace fractals {
using RGB = std::int32_t;
constexpr RGB make_rgb(int r, int g, int b) { return (r << 16) | (g << 8) | b; }

inline int red(RGB i) { return 0xff & (i >> 16); }
inline int green(RGB i) { return 0xff & (i >> 8); }
inline int blue(RGB i) { return 0xff & i; }

} // namespace fractals