#include "complex.hpp"
#include "fractal.hpp"
#include "high_precision_real.hpp"
#include "mandelbrot.hpp"

#include <vector>
#undef NDEBUG
#include <cassert>
#include <iostream>

using namespace fractals;

int main() {
  using namespace mandelbrot;
  using C = std::complex<double>;

  assert((C{0, 0} + C{1, 1} == C{1, 1}));
  assert((square(C{0, 1}) == C{-1, 0}));

  // High precision
  using H1 = high_precision_real<2>;

  H1 h1;
  assert(h1.to_double() == 0);
  H1 h2(12.0);
  assert(h2.to_double() == 12.0);

  h2 = H1{4};
  assert(h2.to_int() == 4);

  H1 h3 = H1{-3};
  assert(h2 * 2 == H1{8});
  assert(h2 / 2 == H1{2});

  assert(H1{1} - H1{2} == H1{-1});
  assert(H1{2} - H1{1} == H1{1});
  assert(H1{2} - H1{-1} == H1{3});
  assert(H1{2} - H1{-3} == H1{5});

  assert(H1{-2} - H1{-3} == H1{1});
  assert(H1{-3} - H1{-2} == H1{-1});

  // Serialisation of high precision
  {
    using H10 = high_precision_real<10>;
    using H2 = high_precision_real<2>;

    std::cout << H2{1.23} << std::endl;

    auto half = H10{0.5};

    for (int i = 0; i < 10; i++) {
      std::cout << half << std::endl;
      half = half * half;
    }

    auto t = H10{1} / 10;
    std::cout << std::hex << t.fraction[1] << std::endl;

    auto x = H10{0};
    x.fraction[1] = 0x1999999999999999ull;
    for (int i = 2; i < 10; ++i) {
      x.fraction[i] = 0x9999999999999999ull;
    }

    std::cout << x << std::endl;
    std::istringstream ss("4.3200000001234567890123456789");

    ss >> x;
    std::cout << x << std::endl;

    // Check the round-trip
    for (int i = 0; i < 100; i++) {
      std::stringstream ss1;
      ss1 << x;
      ss1 >> x;
      std::cout << x << std::endl;
    }
  }
  return 0;
}
