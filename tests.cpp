#include "complex.hpp"
#include "fractal.hpp"
#include "high_precision_real.hpp"
#include "mandelbrot.hpp"
#include "orbit.hpp"
#include "rendering_sequence.hpp"

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

class test_rendering_sequence
    : public fractals::buffered_rendering_sequence<double> {
public:
  test_rendering_sequence()
      : fractals::buffered_rendering_sequence<double>(500, 500, 16) {}

  std::atomic<int> points = 0;
  int layers = 0;
  double get_point(int x, int y) override {
    ++points;
    return x + y;
  }

  void layer_complete(int stride) override { ++layers; }
};

int main() {
  using MB2 = mandelbrot_calculation<2>;
  using C = std::complex<double>;

  auto o1 = make_basic_orbit<MB2>(C{0.25, 0.5});
  compare_orbits(o1, o1, 100);

  {
    fractals::rendering_sequence seq1(5, 10, 16);
    int count = 0;
    int x, y, s;
    bool c;
    while (seq1.next(x, y, s, c))
      ++count;
    assert(count == 50);
  }

  {
    test_rendering_sequence s1;
    std::atomic<bool> stop;
    s1.calculate(1, stop);
    assert(s1.points == 500 * 500);
    assert(s1.layers == 5);
  }
  return 0;
}
