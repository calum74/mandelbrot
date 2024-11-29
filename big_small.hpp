#pragma once

namespace mandelbrot {
template <typename C> struct big {
  using value_type = typename C::value_type;
  C c;

  value_type real_part() const { return c.real_part(); }
  value_type imaginary_part() const { return c.imaginary_part(); }
};

template <typename C> struct small {
  using value_type = typename C::value_type;
  C c;
  value_type real_part() const { return c.real_part(); }
  value_type imaginary_part() const { return c.imaginary_part(); }
};

template <typename B, typename S> struct big_small {
  using value_type = typename B::value_type;
  using big_type = B;
  using small_type = S;

  big_type b;
  small_type s;
  explicit big_small(big_type b) : b{b} {}
  explicit big_small(small_type s) : s{s} {}
  big_small(big_type b0, small_type s0) : b{b0}, s{s0} {
    // Check for potential overflow of s
    if (abs(s.real_part()) > 10 || abs(s.imaginary_part()) > 10) {
      b.c = b.c + s.c;
      s.c = {};
    }
  }

  value_type real_part() const { return b.real_part() + s.real_part(); }
  value_type imaginary_part() const {
    return b.imaginary_part() + s.imaginary_part();
  }
};

template <typename R>
using big_small_t = big_small<big<complex<R>>, small<complex<R>>>;

template <typename C> small<C> operator+(small<C> a, small<C> b) {
  return {a.c + b.c};
}

template <typename C>
big_small<big<C>, small<C>> operator+(big<C> a, small<C> b) {
  return {a, b};
}

template <typename C>
big_small<big<C>, small<C>> operator+(big_small<big<C>, small<C>> a,
                                      big_small<big<C>, small<C>> b) {
  return {a.b + b.b, a.s + b.s};
}

// template <typename C> small<C> operator-(small<C> a, small<C> b) {
//   return {a.c - b.c};
// }

template <typename C> small<C> operator*(small<C> a, small<C> b) {
  return {a.c * b.c};
}

template <typename C> small<C> operator*(int a, small<C> b) {
  return {a * b.c};
}

template <typename C> big<C> operator+(big<C> a, big<C> b) {
  return {a.c + b.c};
}

// template <typename C> big<C> operator-(big<C> a, big<C> b) {
//   return {a.c - b.c};
// }

template <typename C> big<C> operator*(big<C> a, big<C> b) {
  return {a.c * b.c};
}

template <typename C> small<C> operator*(big<C> a, small<C> b) {
  return {a.c * b.c};
}

template <typename C> big<C> operator*(int a, big<C> b) { return {a * b.c}; }

template <typename B, typename S> big_small<B, S> square(big_small<B, S> c) {
  return {square(c.b), 2 * c.b * c.s + square(c.s)};
}

template <typename T> big<T> square(big<T> c) { return {square(c.c)}; }

template <typename T> small<T> square(small<T> c) { return {square(c.c)}; }

} // namespace mandelbrot