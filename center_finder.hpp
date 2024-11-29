#pragma once
#include <vector>

namespace fractals {
class Viewport;

// An algorithm to locate potential centers
//
std::vector<std::pair<int, int>> find_centers(/* const */ Viewport &vp,
                                              int window_size);

} // namespace fractals