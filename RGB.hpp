#pragma once
#include <cstdint>

namespace fractals {
using RGB = std::int32_t;
inline RGB make_rgb(int r, int g, int b) { return (r << 16) | (g << 8) | b; }
inline RGB make_rgbx(int r, int g, int b, int x) {
  return (x << 24) | (r << 16) | (g << 8) | b;
}
inline int red(RGB i) { return 0xff & (i >> 16); }
inline int green(RGB i) { return 0xff & (i >> 8); }
inline int blue(RGB i) { return 0xff & i; }
inline int extra(RGB i) { return 0xff & ((std::uint32_t)i >> 24); }
inline RGB with_extra(RGB i, int x) { return i & (0xffffff) | (x << 24); }

} // namespace fractals