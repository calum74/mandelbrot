/*
    An orbit is an infinite sequence of complex numbers modelled as a C++
    forward iterator or random-access iterator.

    The first element of an orbit is 0.

    The current value of the orbit is given by operator*()
    The orbit is advanced using operator++().
    For random-access orbits, the value of an element is given by operator[].

    I did not invent any of these algorithms, I've just wrapped them in C++
    iterators.
*/

#include "complex.hpp"
#include <vector>
#include <atomic>


namespace mandelbrot {

// !! Move these into complex.hpp
template <int Order, typename C, bool is_even = (Order % 2 == 0)>
struct pow_impl;

template <typename C> struct pow_impl<2, C, true> {
  static C eval(const C &c) { return square(c); }
};

template <typename C> struct pow_impl<1, C, false> {
  static C eval(const C &c) { return c; }
};

template <int Order, typename C> struct pow_impl<Order, C, true> {
  static C eval(const C &c) {
    auto r = pow_impl<Order / 2, C>::eval(c);
    return r * r;
  }
};

template <typename C> struct pow_impl<0, C, true> {
  static C eval(const C &c) { return C{1, 0}; }
};

template <int Order, typename C> struct pow_impl<Order, C, false> {
  static C eval(const C &c) { return c * pow_impl<Order - 1, C>::eval(c); }
};

template <int Order, typename C> C pow(const C &c) {
  return pow_impl<Order, C>::eval(c);
}

template <int N, int M> struct choose_impl {
  static const int value =
      choose_impl<N - 1, M - 1>::value + choose_impl<N - 1, M>::value;
};

template <int M> struct choose_impl<0, M> {
  static const int value = 0;
};

template <> struct choose_impl<0, 1> {
  static const int value = 0;
};

template <int N> struct choose_impl<N, N> {
  static const int value = 1;
};

template <int N> struct choose_impl<N, 0> {
  static const int value = 1;
};

template <int N> struct choose_impl<N, 1> {
  static const int value = N;
};

template <int N, int M> constexpr int choose() {
  return choose_impl<N, M>::value;
}

// The complex arithmetic required to calculate a Mandelbrot set
// We consider the generalized Mandelbrot set, including higher orders
// such as the cubic, but exclude non-integer orders.
// TODO: Specialise this for N=2 for speed
template <int N> struct mandelbrot_calculation {

  // The general form of the calculation z -> z^N + c
  template <typename Complex>
  static Complex step(const Complex &z, const Complex &c) {
    return pow<N>(z) + c;
  }

  template <int J> struct calculate_epsilon {
    template <typename Complex>
    static Complex eval(const Complex &z, const Complex &e) {
      return choose<N, J>() * pow<J>(z) * pow<N - J>(e) +
             calculate_epsilon<J + 1>::eval(z, e);
    }
  };

  template <> struct calculate_epsilon<N> {
    template <typename Complex>
    static Complex eval(const Complex &z, const Complex &e) {
      return {0, 0};
    }
  };

  // When performing perturbations (for higher precision), here is the general
  // formula for evaluating the epsilon (dz)
  template <typename Complex>
  static Complex step_epsilon(const Complex &z, const Complex &e,
                              const Complex &d) {
    return d + calculate_epsilon<0>::eval(z, e);
  }

  template <typename Complex>
  static Complex A(const Complex &z, const Complex &A_prev) {
    return choose<N, 1>() * pow<N - 1>(z) * A_prev + Complex{1, 0};
  }

  template <typename Complex>
  static Complex B(const Complex &z, const Complex &A_prev,
                   const Complex &B_prev) {
    return choose<N, 1>() * pow<N - 1>(z) * B_prev +
           choose<N, 2>() * pow<N - 2>(z) * pow<2>(A_prev);
  }

  template <typename Complex>
  static Complex C(const Complex &z, const Complex &A_prev,
                   const Complex &B_prev, const Complex &C_prev) {
    return choose<N, 1>() * pow<N - 1>(z) * C_prev +
           choose<N, 2>() * pow<N - 2>(z) * Complex{2} * A_prev * B_prev +
           (N > 2 ? choose<N, 3>() * pow<N - 3>(z) * pow<3>(A_prev) : 0);
  }

  // TODO: Can generalize this further to calculate the Nth Taylor series term
};

/*
  A basic_orbit iterates the escape sequence manually.
*/
template <typename Complex, typename Calculation> class basic_orbit {
public:
  using value_type = Complex;

  basic_orbit(value_type c = value_type{}) : c{c}, z{} {}

  const value_type &operator*() const { return z; }

  basic_orbit &operator++() {
    z = Calculation::step(z, c);
    return *this;
  }

private:
  value_type c, z;
};

template <typename Calculation, typename Complex>
basic_orbit<Complex, Calculation> make_basic_orbit(const Complex &c) {
  return {c};
}

/*
  An orbit which is converts its values to a different datatype.
  Generally useful when exposing a high-precision orbit into
  a lower-precision form.

  Complex is a lower-precision complex number (e.g. std::complex<double>)
  ReferenceOrbit is a high precision orbit (e.g.
  basic_orbit<std::complex<high_precision_real<>>>)
*/
template <typename Complex, typename ReferenceOrbit> class converted_orbit {
public:
  converted_orbit() = default;
  converted_orbit(const ReferenceOrbit &r) : reference(r) {}

  using value_type = Complex;

  value_type operator*() const {
    return convert_complex<value_type>(*reference);
  }

  converted_orbit &operator++() {
    ++reference;
    return *this;
  }

private:
  ReferenceOrbit reference;
};

/*
  An orbit which is calculated internally to a high precision, but yields a
  low-precision sequence. This is because most algorithms don't need the
  high-precision number.
 */
template <typename LowPrecisionComplex, typename HighPrecisionComplex,
          typename Calculation>
using high_precision_orbit =
    converted_orbit<LowPrecisionComplex,
                    basic_orbit<HighPrecisionComplex, Calculation>>;

template <typename C, typename Ref, typename Calculation> class relative_orbit {
public:
  using value_type = C;
  relative_orbit(Ref ref, value_type delta)
      : reference_orbit(ref), delta{delta}, epsilon{} {}

  value_type operator*() const { return *reference_orbit + epsilon; }

  value_type reference_z() const { return *reference_orbit; }

  relative_orbit &operator++() {
    epsilon = Calculation::step_epsilon(*reference_orbit, epsilon, delta);
    ++reference_orbit;
    return *this;
  }

  value_type get_epsilon() const { return epsilon; }

private:
  Ref reference_orbit;
  value_type delta, epsilon;
};

  template <typename C, typename Ref> class stored_orbit {
  public:
    stored_orbit() : current{0} { push_next(); }
    stored_orbit(Ref o) : orbit(o), current{0} { push_next(); }

    using value_type = C;

    void push_next() {
      values.push_back(convert_complex<C>(*orbit));
      ++orbit;
    }

    value_type operator*() const { return values[current]; }

    const value_type &operator[](int n) {
      while (values.size() < n)
        push_next();
      return values[n];
    }

    stored_orbit &operator++() {
      ++current;
      if (current >= values.size())
        push_next();
      return *this;
    }

  private:
    Ref orbit;
    std::vector<value_type> values;
    int current;
  };

  template <typename Calculation, typename Rel>
  relative_orbit<typename Rel::value_type, Rel, Calculation>
  make_relative_orbit(Rel rel, auto delta) {
    return {rel, delta};
  }

  /*
    An orbit that also computes terms of a Taylor series approximation for
    each iteration.

    z = z_0 + A.delta + B.delta^2 + C.delta^3
  */
  template <typename Complex, typename ReferenceOrbit, typename Calculation>
  class taylor_series_orbit {
  public:
    using value_type = Complex;
    value_type A, B, C;

    taylor_series_orbit() = default;

    taylor_series_orbit(ReferenceOrbit o) : orbit(o) {
      A = {0.0, 0.0};
      B = {0.0, 0.0};
      C = {0.0, 0.0};
    }

    // Convert the high-precision value into a low-precision value
    value_type operator*() const { return convert_complex<value_type>(*orbit); }

    taylor_series_orbit &operator++() {

      // See
      // https://en.wikipedia.org/wiki/Plotting_algorithms_for_the_Mandelbrot_set
      auto A_next = Calculation::A(**this, A);
      auto B_next = Calculation::B(**this, A, B);
      auto C_next = Calculation::C(**this, A, B, C);

      A = A_next;
      B = B_next;
      C = C_next;

      ++orbit;
      return *this;
    }

  private:
    ReferenceOrbit orbit;
  };

  // An orbit that's relative to another reference orbit, so can be computed
  // using a low-precision complex number. ReferenceOrbit must be a
  // random-access orbit (supporting [])
  template <typename Complex, typename ReferenceOrbit, typename Calculation>
  class perturbation_orbit {
  public:
    using value_type = Complex;

    perturbation_orbit(const ReferenceOrbit &ref, Complex delta,
                       int starting_iteration = 0,
                       Complex starting_epsilon = {})
        : n{starting_iteration}, j{n}, delta(delta), epsilon{starting_epsilon},
          reference{ref}, skipped{starting_iteration} {}

    int iteration() const { return n; }

    value_type operator*() const { return reference[j] + epsilon; }

    perturbation_orbit &operator++() {

      epsilon = Calculation::step_epsilon(reference[j], epsilon, delta);
      j++;

      auto z = reference[j] + epsilon;

      if (escaped(reference[j]) || norm(z) < norm(epsilon)) {
        // We have exceeded the bounds of the current orbit
        // We need to reset the current orbit.
        // Thanks to
        // https://github.com/ImaginaFractal/Algorithms/blob/main/Perturbation.cpp
        // https://philthompson.me/2022/Perturbation-Theory-and-the-Mandelbrot-set.html
        epsilon = z;
        j = 0;
      }

      n++;
      return *this;
    }

    // For information purposes only
    int skipped_iterations() const { return skipped; }

  private:
    int n, j;
    Complex delta, epsilon;
    const ReferenceOrbit &reference;
    int skipped;
  };

  /*
    A Taylor series orbit, which is a reference orbit together with the Taylor
    series terms for each iteration. This allows for an efficient
    (random-access) test. Complex is a lower-precision representation, e.g.
    std::complex<double> ReferenceOrbit is a high-precision orbit.
  */
  template <typename Complex, typename ReferenceOrbit, typename Calculation>
  class stored_taylor_series_orbit {
  public:
    using value_type = Complex;

    stored_taylor_series_orbit() = default;
    stored_taylor_series_orbit(ReferenceOrbit r) : reference(r) {}
    stored_taylor_series_orbit(ReferenceOrbit r, int max_iterations,
                               std::atomic<bool> &stop)
        : reference(r) {
      for (int i = 0; !stop && i <= max_iterations && !escaped((*this)[i]); ++i)
        ; // Force evaluation of the reference orbit
    }

    // Gets the corresponding epsilon (dz) for a given delta (dc)
    // at iteration i.
    auto epsilon(int i, Complex delta) {
      while (i >= entries.size())
        get_next();
      return entries[i].epsilon(delta);
    }

    auto epsilon(int i, Complex delta) const {
      return entries.at(i).epsilon(delta);
    }

    // Gets the value of the reference orbit at iteration i
    Complex operator[](int i) {
      while (i >= entries.size())
        get_next();
      return entries[i].z;
    }

    Complex operator[](int i) const { return entries.at(i).z; }

  private:
    // Finds the number of iterations it's safe to skip
    // because we judge that we are sufficiently close to the reference orbit
    std::pair<int, Complex> find_iterations_to_skip(Complex delta,
                                                    int max) const {
      int min = 0;

      // Quoting Superfractalthing Maths by K.I. Martin,
      // The approximation should be good as long as the δ^3
      // term has a magnitude significantly smaller then the δ^2 term.
      //
      // This means that we'll use the approximation until the terms run out,
      // then we transition to iteration.

      Complex epsilon = {};
      while (max - min > 8) {
        int mid = (max + min) / 2;
        auto e = this->epsilon(mid, delta);
        if (e.second && !escaped((*this)[mid])) {
          min = mid;
          epsilon = e.first;
        } else
          max = mid;
      }

      return {min, epsilon};
    }

  private:
    struct Entry {
      Complex z, A, B, C;

      // Returns the epsilon, and whether the epsilon is "accurate"
      std::pair<Complex, bool> epsilon(Complex delta) const {
        auto delta_2 = delta * delta;
        auto t2 = B * delta_2;
        auto t3 = C * delta * delta_2;
        // TODO: Can we precalculate if this entry is accurate?
        return std::make_pair(A * delta + t2 + t3, norm(t2) > 10000 * norm(t3));
      }
    };

    taylor_series_orbit<Complex, ReferenceOrbit, Calculation> reference;
    std::vector<Entry> entries;

    void get_next() {
      entries.push_back({*reference, reference.A, reference.B, reference.C});
      ++reference;
    }

  public:
    using relative_orbit =
        perturbation_orbit<value_type, stored_taylor_series_orbit, Calculation>;

    relative_orbit make_relative_orbit(Complex delta, int limit) const {
      auto s = find_iterations_to_skip(delta, entries.size());
      return {*this, delta, s.first, s.second};
    }
  };

} // namespace mandelbrot
