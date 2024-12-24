#include "high_precision_real.hpp"
#undef NDEBUG
#include <cassert>
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
  if (a + b - b != a) {
    auto dbg1 = a + b;
    auto dbg2 = dbg1 - b;
  }
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
    // TODO: This fails when repeated 10000 times
    assert_eq(a * inverse(a), one);
    assert_eq(inverse(a) * a, one);

    // TODO: Division

    a.fraction[0] = 0;
    test_arithmetic(a, b);
    test_arithmetic(b, a);
    a.fraction[1] = 0;
    test_arithmetic(a, b);
    test_arithmetic(b, a);
  }
}

int main() {

  // Initializers
  assert(hp<2>{}.to_double() == 0);

  test<2>(100);
  test<3>(100);
  test<4>(100);
  test<10>(100);

  // To string
  assert(to_string(hp<2>{}) == "0.000000");
}
