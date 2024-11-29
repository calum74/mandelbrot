#include "center_finder.hpp"
#include "RGB.hpp"
#include "Viewport.hpp"

#include <algorithm>
#include <cassert>
#include <iostream>

namespace {
using namespace fractals;
struct rgb_sum {
  long sum = 0;

  rgb_sum &operator+=(RGB p) {
    sum += red(p);
    sum += green(p);
    sum += blue(p);
    return *this;
  }
  rgb_sum &operator-=(RGB p) {
    sum -= red(p);
    sum -= green(p);
    sum -= blue(p);
    return *this;
  }
};
struct pixel_data {
  long votes = 0;
  rgb_sum sum[4];
};

enum direction { LR, RL, UD, DU };
}; // namespace

std::vector<std::pair<int, int>> fractals::find_centers(Viewport &vp,
                                                        int window_size) {

  std::vector<std::pair<int, int>> result;

  return result;

  /*
  Algorithm: Locate an image on the left hand side, and see if we can find the
  same image in the image somewhere.

  Algorithm: For every point, get a measure of its "symmetry" to a certain size.

  */

  auto difference = [](RGB a, RGB b) {
    auto d1 = red(a) - red(b);
    auto d2 = green(a) - green(b);
    auto d3 = blue(a) - blue(b);
    return d1 * d1 + d2 * d2 + d3 * d3;
  };

  auto compute_symmetry = [&](int x, int y, double bailout) {
    if (x >= window_size && x < vp.width - window_size && y >= window_size &&
        y <= vp.height - window_size) {
      double total = 0;
      for (int j = -window_size; total < bailout && j < window_size; ++j)
        for (int i = 0; total < bailout && i < window_size; ++i) {
          auto p1 = vp(x + i, y - j);
          auto p2 = vp(x - i, y + j);
          total += difference(p1, p2);
        }
      return total;
    } else {
      return 1e50;
    }
  };

  double best = 1e50;
  int best_x = 0, best_y = 0;

  for (int j = window_size; j < vp.height - window_size; ++j)
    for (int i = window_size; i < vp.width - window_size; ++i) {
      auto d = compute_symmetry(i, j, best);
      if (d < best) {
        best_x = i;
        best_y = j;
        best = d;
      }
    }
  vp(best_x, best_y) = make_rgb(255, 0, 0);

  std::cout << "Best = " << best_x << "," << best_y << std::endl;
}