#include "complex_number.hpp"
#include "fractal.hpp"
#include "high_exponent_real.hpp"
#include "high_precision_real.hpp"
#include "mandelbrot.hpp"
#include "orbit.hpp"
#include "orbit_manager.hpp"
#include "orbit_tree.hpp"
#include "rendering_sequence.hpp"
#include "view_parameters.hpp"

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
    if (escaped(*o1))
      return;
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

void test_sequence(int w, int h) {
  std::vector<bool> visited(w * h);
  for (int i = 0; i < w * h; ++i) {
    multi_resolution_sequence p(w, h, i);
    visited.at(p.x + p.y * w) = true;
  }
  for (bool b : visited) {
    if (!b) {
      for (int i = 0; i < w * h; ++i) {
        multi_resolution_sequence p(w, h, i);
        std::cout << i << " (" << p.x << "," << p.y << ")\n";
      }
    }
    assert(b);
  }
}

int main() {
  {
    // Basic rendering sequence
    test_sequence(3, 4);
    test_sequence(4, 4);
    for (int w = 0; w < 9; ++w)
      for (int h = 0; h < 9; ++h)
        test_sequence(w, h);
  }

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

  {
    fractals::view_coords c1, c2;
    c1.x = 1.5;
    c1.y = -0.5;
    c1.r = 0.25;
    c1.max_iterations = 500;
    std::stringstream ss1, ss2;
    ss1 << c1;
    ss1 >> c2;
    ss2 << c2;
    assert(ss2.str() == "(1.5000,-0.5000,0.2499,500)");
  }

  {
    fractals::view_parameters p1{.coords = {1.5, -0.5, 0.5, 500},
                                 .fractal_name = "abc def",
                                 .colour_seed = 99,
                                 .colour_gradient = 0.001},
        p2;
    std::stringstream ss1, ss2;
    ss1 << p1;
    ss1 >> p2;
    ss2 << p2;
    assert(ss1.str() == ss2.str());
    assert(p1.fractal_name == p2.fractal_name);
    assert(p1.colour_seed == p2.colour_seed);
  }

  // Test Mandelbrot deltas

  {
    using Calc = mandelbrot::mandelbrot_calculation<2>;
    std::complex<double> c = {-0.53235, -0.60034}, A = 0, B = 0, C = 0, D = 0,
                         z = 0;

    c = {-0.2238286049999855391525793, -1.1167864957492792902234038};
    constexpr int N = 10;
    std::array<std::complex<double>, N> V, V2, zero;

    for (int i = 0; i < 500 && !escaped(z); i++) {
      auto A2 = Calc::A(z, A);
      auto B2 = Calc::B(z, A, B);
      auto C2 = Calc::C(z, A, B, C);
      auto D2 = Calc::D(z, A, B, C, D);
      A = A2;
      B = B2;
      C = C2;
      D = D2;
      V = Calc::delta_terms(z, V);
      z = Calc::step(z, c);

      /*
      std::cout << z << "   " << A << B << C << "   ";
      for (int j = 0; j < N; ++j)
        std::cout << V[j];
      std::cout << std::endl;
      */

      auto approx_eq = [](std::complex<double> a, std::complex<double> b) {
        return norm(a) == 0 && norm(b) == 0 || norm(a / b - 1.0) < 0.0001;
      };
      assert(approx_eq(A, V[0]));
      assert(approx_eq(B, V[1]));
      assert(approx_eq(C, V[2]));
      assert(approx_eq(D, V[3]));
    }
  }

  {
    mandelbrot::basic_orbit<std::complex<double>,
                            mandelbrot::mandelbrot_calculation<2>>
        orbit1({0.5, 0.5});

    int iterations = 0;
    while (!mandelbrot::escaped(*orbit1)) {
      ++iterations;
      ++orbit1;
    }

    assert(iterations == 5);
  }

  {
    using R1 = fractals::high_exponent_real<double, int>;
    using R2 = double;
    using R3 = fractals::high_precision_real<3>;

    mandelbrot::basic_orbit<std::complex<R1>,
                            mandelbrot::mandelbrot_calculation<2>>
        orbit1({0.5, 0.5});
    mandelbrot::basic_orbit<std::complex<R2>,
                            mandelbrot::mandelbrot_calculation<2>>
        orbit2({0.5, 0.5});
    mandelbrot::basic_orbit<std::complex<R3>,
                            mandelbrot::mandelbrot_calculation<2>>
        orbit3({0.5, 0.5});

    int iterations = 0;
    while (!mandelbrot::escaped(*orbit1)) {
      ++iterations;
      ++orbit1;
      ++orbit2;
      ++orbit3;
    }

    std::cout << *orbit1 << *orbit2 << *orbit3;

    assert(iterations == 5);

    auto reference_orbit =
        mandelbrot::make_basic_orbit<mandelbrot::mandelbrot_calculation<2>>(
            std::complex<R1>{0.49, 0.49});

    std::atomic<bool> stop;
    auto stored_orbit = mandelbrot::make_stored_orbit<std::complex<R1>>(
        reference_orbit, 100, stop);

    mandelbrot::stored_taylor_series_orbit<
        std::complex<R1>, std::complex<R1>, std::complex<R1>,
        mandelbrot::basic_orbit<std::complex<R1>,
                                mandelbrot::mandelbrot_calculation<2>>,
        3, 20, 100>
        taylor_series{reference_orbit, 100, stop};

    auto relative =
        taylor_series.make_relative_orbit({0.01, 0.01}, 100, iterations);
    iterations = relative.iteration();

    // Probably this will change
    // assert(iterations == 3);

    // Iterate the relative orbit as before, using perturbation.
    while (!mandelbrot::escaped(*relative) && iterations < 100) {
      ++iterations;
      ++relative;
    }

    assert(iterations == 5);
  }

  {

    mandelbrot::orbit_manager<
        std::complex<double>,
        std::complex<fractals::high_exponent_real<double>>,
        std::complex<fractals::high_exponent_real<double>>, 4, 20, 100,
        mandelbrot::basic_orbit<std::complex<double>,
                                mandelbrot::mandelbrot_calculation<2>>>
        manager;

    std::atomic<bool> stop;
    manager.new_view({0, 0}, {1, 1}, 4, std::complex{0.5, 0.5}, 100, stop);

    auto orbit1 = manager.lookup({0.1, 0.1}, 100);

    compare_orbits(
        orbit1,
        mandelbrot::make_basic_orbit<mandelbrot::mandelbrot_calculation<2>>(
            std::complex{0.6, 0.6}),
        100);

    auto orbit2 = manager.lookup({0.1, 0.1}, 100);

    compare_orbits(
        orbit2,
        mandelbrot::make_basic_orbit<mandelbrot::mandelbrot_calculation<2>>(
            std::complex{0.6, 0.6}),
        100);

    auto orbit3 = manager.lookup({-1.0, -1.0}, 100);

    // Let's move the reference orbit
    manager.new_view({-0.1, -0.2}, {1, 1}, 100, std::complex{0.4, 0.3}, 100,
                     stop);

    auto orbit4 = manager.lookup({0.2, 0.3}, 100);
    compare_orbits(
        orbit4,
        mandelbrot::make_basic_orbit<mandelbrot::mandelbrot_calculation<2>>(
            std::complex{0.6, 0.6}),
        100);

    // Compare the new orbits
    auto orbit5 = manager.lookup({0.2, 0.3}, 100);
    compare_orbits(
        orbit5,
        mandelbrot::make_basic_orbit<mandelbrot::mandelbrot_calculation<2>>(
            std::complex{0.6, 0.6}),
        100);
  }

  // Orbit-trees
  {
    using ref_type =
        mandelbrot::basic_orbit<std::complex<double>,
                                mandelbrot::mandelbrot_calculation<2>>;
    using stored_type = stored_orbit<std::complex<double>, ref_type>;
    using tree_type =
        mandelbrot::orbit_branch<std::complex<double>, std::complex<double>,
                                 std::complex<double>, stored_type, 4, 10, 10>;
    ref_type r0({-1.248193761, 0.089224601});
    std::atomic<bool> stop;
    stored_type stored(r0, 500, stop);

    std::cout << stored.size() << " iterations in the reference orbit\n";

    std::complex<double> radius = {0.000001275, 0.000001275};
    auto root = std::make_shared<tree_type>(stored, radius, 500, stop);
    std::cout << root->size() << " iterations in the root branch\n";

    auto branch1 = std::make_shared<tree_type>(
        root, mandelbrot::topleft(radius), 500, stop);
    auto branch2 = std::make_shared<tree_type>(
        root, mandelbrot::topright(radius), 500, stop);
    auto branch3 = std::make_shared<tree_type>(
        root, mandelbrot::bottomleft(radius), 500, stop);
    auto branch4 = std::make_shared<tree_type>(
        root, mandelbrot::bottomright(radius), 500, stop);
    std::cout << branch1->size() << " iterations in branch1\n";
    std::cout << branch2->size() << " iterations in branch2\n";
    std::cout << branch3->size() << " iterations in branch3\n";
    std::cout << branch4->size() << " iterations in branch4\n";

    auto x = branch2->get_escape_iterations({0, 0}, 500);
    x = branch2->get_escape_iterations(radius * 0.1, 500);
    std::cout << x << std::endl;

    // orbit_tree<std::complex<double>, std::complex<double>,
    // std::complex < double >> X;
  }

  return 0;
}
