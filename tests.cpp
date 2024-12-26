#include "complex.hpp"
#include "fractal.hpp"
#include "high_precision_real.hpp"
#include "mandelbrot.hpp"
#include "orbit.hpp"

#include <iomanip>
#include <vector>

#undef NDEBUG
#include <cassert>
#include <iostream>

using namespace fractals;
using namespace mandelbrot;

template <typename Orbit1, typename Orbit2>
void compare_orbits(Orbit1 o1, Orbit2 o2, int n) {
  for (int i = 0; i < n; ++i) {
    auto d = norm(*o1 - *o2);
    assert(d < 0.0001);
    ++o1;
    ++o2;
  }
}

int main() {
  using MB2 = mandelbrot_calculation<2>;
  using C = std::complex<double>;

  auto o1 = make_basic_orbit<MB2>(C{0.25, 0.5});
  compare_orbits(o1, o1, 100);
  return 0;
}
