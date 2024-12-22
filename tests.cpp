#include "complex.hpp"
#include "fractal.hpp"
#include "high_precision_real.hpp"
#include "mandelbrot.hpp"

#include <iomanip>
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
  using H2 = high_precision_real<2>;

  {
    H2 h1 = 4;
    auto h3dbg = h1 * 2;
    assert(h1 * 2 == H2{8});
  }

  {
    H2 a = 0.5, b = 1.5;
    std::cout << a * b << ' '; // < (D2 * H2{32.0 / 17.0}) << " ";
  }

  H2 h1;
  assert(h1.to_double() == 0);
  H2 h2(12.0);
  assert(h2.to_double() == 12.0);

  h2 = H2{4};
  assert(h2.to_int() == 4);

  H2 h3 = H2{-3};
  auto h3dbg = h2 * 2;
  assert(h2 * 2 == H2{8});
  assert(h2 / 2 == H2{2});

  assert(H2{1} - H2{2} == H2{-1});
  assert(H2{2} - H2{1} == H2{1});
  assert(H2{2} - H2{-1} == H2{3});
  assert(H2{2} - H2{-3} == H2{5});

  assert(H2{-2} - H2{-3} == H2{1});
  assert(H2{-3} - H2{-2} == H2{-1});

  // Serialisation of high precision
  {
    using H220 = high_precision_real<10>;

    std::cout << H2{1.23} << std::endl;

    auto half = H220{0.5};

    for (int i = 0; i < 10; i++) {
      std::cout << half << std::endl;
      half = half * half;
    }

    auto t = H2{1} / 10;
    std::cout << std::hex << t.fraction[1] << std::endl;

    auto x = H2{0};
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
      ss1 << std::setprecision(100) << x;
      ss1 >> x;
      std::cout << std::dec << x << std::endl;
    }
  }

  // Subtraction
  {
    H2 a = 2.8;
    H2 b = 0.9;
    std::cout << "a-b=" << (a - b) << std::endl;
  }

  {
    H2 x{1.5};
    std::cout << (x << 1) << " " << (x >> 1) << " ";
    std::cout << inverse(H2{1 / 1.97254});
  }
  return 0;
}
