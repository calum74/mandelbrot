#include "high_precision_real.hpp"
#undef NDEBUG
#include <cassert>
#include <iomanip>
#include <random>
#include <sstream>

template <int N> using hp = fractals::high_precision_real<N>;

template <int N>
void assert_eq(const hp<N> &a, const hp<N> &b,
               std::uint64_t tolerance = 0xffff) {
  auto d = a - b;

  for (int i = 0; i < N - 1; i++) {
    assert(d.fraction[i] == 0);
  }
  assert(d.fraction[N - 1] <= tolerance);
}

template <int N> std::string to_string(const hp<N> &n) {
  std::stringstream ss;
  ss << n;
  return ss.str();
}

template <int N> hp<N> make_random(std::default_random_engine &g) {
  hp<N> result;
  std::uniform_int_distribution<std::uint64_t> j;
  for (int i = 1; i < N; ++i)
    result.fraction[i] = j(g);
  j = std::uniform_int_distribution<std::uint64_t>(0, 1024);
  result.fraction[0] = std::uniform_int_distribution<std::uint64_t>(0, 1024)(g);
  result.negative = j(g) & 1;
  return result;
}

template <int N> void test_arithmetic(const hp<N> &a, const hp<N> &b) {
  const hp<N> zero, one{1};

  assert(a + b == b + a);
  assert(a - a == zero);
  assert(((a + a) >> 1) == a);
  assert(a - b == -(b - a));
  assert(a + b - b == a);

  assert(a == -(-a));
  assert(a + zero == a);
  assert(zero - a == -a);
  assert(a - zero == a);
  assert(zero - a == -a);

  // Shift
  assert((a << 1) >> 1 == a);
  assert((a << 3) >> 3 == a);

  // Multiplication
  assert(a * b == b * a);
  assert(-a * b == -b * a);
  assert(a * b == -b * -a);

  assert(one * a == a);
  assert(a * one == a);
  assert(zero * a == zero);
  assert(a * zero == zero);
}

template <int N> void test(int iterations) {
  std::default_random_engine g;

  for (int i = 0; i < iterations; i++) {
    auto a = make_random<N>(g);
    auto b = make_random<N>(g);
    const hp<N> zero, one{1};

    test_arithmetic(a, b);

    // Inverse
    if (a >= 1)
    {
        // Fails for very small a
        assert_eq(a * inverse(a), one);
        assert_eq(inverse(a) * a, one);
    }
    // TODO: Division

    a.fraction[0] = 0;
    test_arithmetic(a, b);
    test_arithmetic(b, a);
    a.fraction[1] = 0;
    test_arithmetic(a, b);
    test_arithmetic(b, a);
  }
}

void random_bits() {
  using namespace fractals;
  // High precision
  using H2 = hp<2>;

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

  using H3 = high_precision_real<3>;
  using H4 = high_precision_real<4>;

  {
    std::cout << inverse(H2{0.506960}) << ' ';

    std::cout << inverse(H4{0.506960}) << ' ';
    H3 tmp;
    tmp.fraction[1] = (std::uint64_t)(-1);
    tmp.fraction[2] = (std::uint64_t)(-1);
    auto failure = H3{1} - tmp;

    assert_eq(failure, H3{0});
    std::cout << failure << ' ';
    // assert(std::tostring(failure) == 0.000000);
  }

  // To string
  assert(to_string(hp<2>{}) == "0.000000");

  // Inverse regression test case
  {
    fractals::high_precision_real<4> n;
    n.fraction[0] = 3;
    n.fraction[1] = 1895646463175238861ull;
    n.fraction[2] = 18205977988746887042ull;
    n.fraction[3] = 6985711868819929312ull;

    auto i = inverse(n);
    auto j1 = n * i;
    auto j2 = i * n;
    assert(j1 == j2);
    fractals::high_precision_real<4> one{1};
    std::cout << std::setprecision(80) << "j1 = " << j1 << std::endl;
    assert_eq(j1, one, 0xffff);
  }

  {
    hp<4> tenth, ten{10}, one{1};
    make_tenth(tenth);
    std::cout << std::dec << std::setprecision(1000) << tenth << std::endl;
    std::cout << (ten * tenth) << std::endl;
    assert_eq(ten * tenth, one);

    // Check string precision
    auto s1 = "1.2589232485349380378421";
    std::stringstream ss1(s1), ss2;
    hp<4> n1, n2;
    ss1 >> n1;
    ss2 << std::setprecision(10000) << n1;
    ss2 >> n2;
    auto debug = n1 - n2;
    assert_eq(n1, n2);
  }
}

int main() {

  // Initializers
  assert(hp<2>{}.to_double() == 0);

  int n = 1000;
  test<2>(n);
  test<3>(n);
  test<4>(n);
  test<10>(n);

  random_bits();
}
