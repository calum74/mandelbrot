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

template <int Order, typename C, bool is_even = (Order % 2 == 0)>
struct pow_impl;

template <typename C> struct pow_impl<2, C, true> {
  static C eval(const C &c) { return c * c; }
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

template <int Order, typename C> struct pow_impl<Order, C, false> {
  static C eval(const C &c) { return c * pow_impl<Order - 1, C>::eval(c); }
};

template <int Order, typename C> C pow(const C &c) {
  return pow_impl<Order, C>::eval(c);
}

template <typename C, int Order> struct algorithm {
  static C step(const C &z, const C &c) { return pow<Order>(z) + c; }

  static C epsilon_step(const C &z, const C &e) { return {}; }
};

/*
  A basic_orbit iterates the escape sequence manually.
*/
template <typename Complex> class basic_orbit {
public:
  using value_type = Complex;

  basic_orbit(value_type c = value_type{}) : c{c}, z{} {}

  const value_type &operator*() const { return z; }

  basic_orbit &operator++() {
    z = step(z, c);
    return *this;
  }

private:
  value_type c, z;
};

template <typename C> basic_orbit<C> make_basic_orbit(const C &c) {
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
template <typename LowPrecisionComplex, typename HighPrecisionComplex>
using high_precision_orbit =
    converted_orbit<LowPrecisionComplex, basic_orbit<HighPrecisionComplex>>;

template <typename C, typename Ref> class relative_orbit {
public:
  using value_type = C;
  relative_orbit(Ref ref, value_type delta)
      : reference_orbit(ref), delta{delta}, epsilon{} {}

  value_type operator*() const { return *reference_orbit + epsilon; }

  value_type reference_z() const { return *reference_orbit; }

  relative_orbit &operator++() {
    epsilon = 2 * *reference_orbit * epsilon + epsilon * epsilon + delta;
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

  template <typename Rel>
  relative_orbit<typename Rel::value_type, Rel>
  make_relative_orbit(Rel rel, auto delta) {
    return {rel, delta};
  }

  /*
    An orbit that also computes terms of a Taylor series approximation for
    each iteration.

    z = z_0 + A.delta + B.delta^2 + C.delta^3
  */
  template <typename Complex, typename ReferenceOrbit>
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
      auto z_2 = 2 * **this;
      auto A_next = z_2 * A + Complex{1.0, 0.0};
      auto B_next = z_2 * B + A * A;
      auto C_next = z_2 * C + 2 * A * B;
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
  template <typename Complex, typename ReferenceOrbit>
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

      epsilon = 2 * reference[j] * epsilon + epsilon * epsilon + delta;
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
  template <typename Complex, typename ReferenceOrbit>
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

    taylor_series_orbit<Complex, ReferenceOrbit> reference;
    std::vector<Entry> entries;

    void get_next() {
      entries.push_back({*reference, reference.A, reference.B, reference.C});
      ++reference;
    }

  public:
    using relative_orbit =
        perturbation_orbit<value_type, stored_taylor_series_orbit>;

    relative_orbit make_relative_orbit(Complex delta, int limit) const {
      auto s = find_iterations_to_skip(delta, entries.size());
      return {*this, delta, s.first, s.second};
    }
  };

} // namespace mandelbrot
