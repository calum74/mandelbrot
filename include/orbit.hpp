/*
    An orbit is an infinite sequence of complex numbers modelled as a C++
    forward iterator or random-access iterator.

    The first element of an orbit is 0.

    The current value of the orbit is given by operator*()
    The orbit is advanced using operator++().
    For random-access orbits, the value of an element is given by operator[].
*/

#include "mandelbrot_calculation.hpp"
#include <atomic>
#include <vector>

namespace mandelbrot {

/*
  A basic_orbit iterates the escape sequence manually.
  Often this is used as a "reference orbit" for more advanced algorithms.

  Complex - the type of the complex number
  Calculation - how to calculate the sequence
*/
template <typename Complex, typename Calculation> class basic_orbit {
public:
  using value_type = Complex;
  using calculation = Calculation;

  basic_orbit(value_type c = value_type{}) : c{c}, z{} {}

  const value_type &operator*() const { return z; }

  basic_orbit &operator++() {
    z = calculation::step(z, c);
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
  using calculation = typename ReferenceOrbit::calculation;

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
  low-precision sequence. Most algorithms don't need the
  high-precision number.
 */
template <typename LowPrecisionComplex, typename HighPrecisionComplex,
          typename Calculation>
using high_precision_orbit =
    converted_orbit<LowPrecisionComplex,
                    basic_orbit<HighPrecisionComplex, Calculation>>;

template <typename C, typename Ref> class relative_orbit {
public:
  using value_type = C;
  using calculation = typename Ref::calculation;

  relative_orbit(Ref ref, value_type delta)
      : reference_orbit(ref), delta{delta}, epsilon{} {}

  value_type operator*() const { return *reference_orbit + epsilon; }

  value_type reference_z() const { return *reference_orbit; }

  relative_orbit &operator++() {
    epsilon = calculation::step_epsilon(*reference_orbit, epsilon, delta);
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
  stored_orbit(Ref o) : orbit(o) {}

  using value_type = C;
  using calculation = typename Ref::calculation;

  const value_type &operator[](int n) {
    while (values.size() <= n)
      push_next();
    return values[n];
  }

  struct reference_orbit {
  public:
    using value_type = C;
    using calculation = typename Ref::calculation;
    reference_orbit(stored_orbit &ref) : ref{ref}, current{0} {}

    const value_type &operator*() const { return ref[current]; }

    reference_orbit &operator++() {
      ++current;
      return *this;
    }

  private:
    int current;

    stored_orbit &ref;
  };

  reference_orbit make_reference() { return {*this}; }

private:
  void push_next() {
    values.push_back(convert_complex<C>(*orbit));
    ++orbit;
  }

  Ref orbit;
  std::vector<value_type> values;
};

template <typename C, typename Ref>
stored_orbit<C, Ref> make_stored_orbit(const Ref &r) {
  return {r};
}

template <typename Rel>
relative_orbit<typename Rel::value_type, Rel> make_relative_orbit(Rel rel,
                                                                  auto delta) {
  return {rel, delta};
}

/*
  An orbit that also computes terms of a Taylor series approximation for
  each iteration.

  z = z_0 + A.delta + B.delta^2 + C.delta^3
*/
template <typename Complex, typename ReferenceOrbit, int Terms>
class taylor_series_orbit {
public:
  using value_type = Complex;
  using calculation = typename ReferenceOrbit::calculation;

  std::array<Complex, Terms> terms;

  taylor_series_orbit() = default;

  taylor_series_orbit(ReferenceOrbit o) : orbit(o) {}

  // Convert the high-precision value into a low-precision value
  value_type operator*() const { return convert_complex<value_type>(*orbit); }

  taylor_series_orbit &operator++() {

    // See
    // https://en.wikipedia.org/wiki/Plotting_algorithms_for_the_Mandelbrot_set
    terms = calculation::delta_terms(**this, terms);
    ++orbit;
    return *this;
  }

private:
  ReferenceOrbit orbit;
};

// An orbit that's relative to another reference orbit, so can be computed
// using a low-precision complex number. ReferenceOrbit must be a
// random-access orbit (supporting [])
template <typename Complex, typename HighExponentComplex,
          typename ReferenceOrbit>
class perturbation_orbit {
public:
  using value_type = Complex;
  using calculation = typename ReferenceOrbit::calculation;
  using delta_type = HighExponentComplex;
  using epsilon_type = delta_type;

  perturbation_orbit(const ReferenceOrbit &ref, delta_type delta,
                     int starting_iteration = 0,
                     epsilon_type starting_epsilon = {})
      : n{starting_iteration}, j{n}, delta(delta), epsilon{starting_epsilon},
        reference{ref} {}

  int iteration() const { return n; }

  value_type operator*() const {
    return reference[j] + convert_complex<value_type>(epsilon);
  }

  perturbation_orbit &operator++() {

    epsilon = calculation::step_epsilon(reference[j], epsilon, delta);
    j++;

    auto z = reference[j] + epsilon;

    if (escaped(reference[j]) || fractals::norm(z) < fractals::norm(epsilon)) {
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

private:
  int n, j;
  delta_type delta;
  epsilon_type epsilon;
  const ReferenceOrbit &reference;
};

/*
  A Taylor series orbit, which is a reference orbit together with the Taylor
  series terms for each iteration. This allows for an efficient
  (random-access) test.

  Complex is a lower-precision representation, e.g.
  `std::complex<double>`

  HighExponentComplex is used for the Taylor series terms, whose exponents can
  often exceed a standard `double` data type.

  ReferenceOrbit is the high-precision orbit.
*/
template <typename Complex, typename HighExponentComplex,
          typename ReferenceOrbit, int Terms, int Precision>
class stored_taylor_series_orbit {
public:
  using value_type = Complex;
  using term_type = HighExponentComplex;
  using delta_type = HighExponentComplex;
  using epsilon_type = HighExponentComplex;
  using calculation = typename ReferenceOrbit::calculation;

  stored_taylor_series_orbit() = default;
  stored_taylor_series_orbit(ReferenceOrbit r) : reference(r) {}
  stored_taylor_series_orbit(ReferenceOrbit r, int max_iterations,
                             std::atomic<bool> &stop)
      : reference(r) {
    for (int i = 0; !stop && i <= max_iterations && !escaped((*this)[i]); ++i)
      ; // Force evaluation of the reference orbit
    // debug_terms();
  }

  void debug_terms() const {
    // Debug - look at the terms
    int max_term[Terms] = {};
    using R = HighExponentComplex::value_type;
    R max[Terms] = {}, min[Terms] = {};
    for (int i = 0; i < entries.size(); ++i) {
      //
      bool failed = false;
      for (int t = 0; t < Terms; ++t) {
        if ((!isfinite(fractals::real_part(entries[i].terms[t])) ||
             !isfinite(fractals::imag_part(entries[i].terms[t])))) {
          if (!max_term[t]) {
            max_term[t] = i;
          }
        } else {
          // It's finite
          auto nt = fractals::norm(entries[i].terms[t]);
          if (max[t] == R(0) || nt > max[t])
            max[t] = nt;
          if (min[t] == R(0) || nt < min[t])
            min[t] = nt;
        }
      }
    }
    for (int t = 0; t < Terms; ++t) {
      if (max_term[t]) {
        std::cout << "Maximum term for " << t << " = " << max_term[t];
        for (int u = 0; u < Terms; ++u)
          std::cout << " " << entries[max_term[t]].terms[u];
        std::cout << std::endl;
      }
      std::cout << "Max = " << max[t] << ", min = " << min[t] << std::endl;
    }
  }

  // Gets the corresponding epsilon (dz) for a given delta (dc)
  // at iteration i.
  auto epsilon(int i, delta_type delta) {
    while (i >= entries.size())
      get_next();
    return entries[i].epsilon(delta);
  }

  auto epsilon(int i, delta_type delta) const {
    return entries.at(i).epsilon(delta);
  }

  // Gets the value of the reference orbit at iteration i
  value_type operator[](int i) {
    while (i >= entries.size())
      get_next();
    return entries[i].z;
  }

  value_type operator[](int i) const { return entries.at(i).z; }

private:
  // Finds the number of iterations it's safe to skip
  // because we judge that we are sufficiently close to the reference orbit
  value_type find_iterations_to_skip(delta_type delta, int max,
                                     int &iterations_skipped) const {
    int min = 0;

    // Step 1: Establish the max and min
    Complex epsilon = {};
    int window_size = 4;

    int skipped = iterations_skipped;

    if (skipped > max)
      skipped = max;

    auto e = this->epsilon(skipped, delta);
    if (e.second && !escaped((*this)[skipped])) {
      // Seek upwards
      min = skipped;
      epsilon = convert_complex<Complex>(e.first);
      for (int mid = iterations_skipped + window_size; mid < max;
           mid += (window_size *= 2)) {
        e = this->epsilon(mid, delta);
        if (e.second && !escaped((*this)[mid])) {
          min = mid;
          epsilon = convert_complex<Complex>(e.first);
        } else {
          max = mid;
          break;
        }
      }
    } else {
      // Seek downwards
      max = skipped;
      for (int mid = iterations_skipped - window_size; mid >= 0;
           mid -= (window_size *= 2)) {
        e = this->epsilon(mid, delta);
        if (e.second && !escaped((*this)[mid])) {
          min = mid;
          epsilon = convert_complex<Complex>(e.first);
          break;
        } else {
          max = mid;
        }
      }
    }

    // Step 2:
    while (max - min > 4) {
      int mid = (max + min) / 2;
      auto e = this->epsilon(mid, delta);
      if (e.second && !escaped((*this)[mid])) {
        min = mid;
        epsilon = convert_complex<Complex>(e.first);
      } else
        max = mid;
    }
    iterations_skipped = min;

    return epsilon;
  }

private:
  // An `Entry` is an orbital value (`z`) and the Taylor series terms (`terms`)
  // for an iteration.
  struct Entry {

    using norm_type = typename term_type::value_type;

    Entry(value_type z, const std::array<term_type, Terms> &terms)
        : z(z), terms(terms) {
      // We'll look at the terms in the series to figure out what the maximum
      // size of delta is for this term before imprecision sets in.
      // Each term has a "norm" giving an indication of its size.
      // We need to make sure that each term is sufficiently "small"
      // relative to the previous term.

      auto prev_norm = fractals::norm(terms[0]);
      for (int i = 1; i < Terms - 1; i++) {
        auto n = fractals::norm(terms[i]);
        auto nr = prev_norm / (norm_type(10) * n);
        prev_norm = n;
        if (i == 1 || max_delta_norm > nr)
          max_delta_norm = nr;
      }
      auto n = fractals::norm(terms[Terms - 1]);
      auto nr = prev_norm / (norm_type(100) * n);

      if (max_delta_norm > nr)
        max_delta_norm = nr;
    }
    value_type z;
    std::array<term_type, Terms> terms;

    // The maximum delta for this term
    // If we use a term with higher delta than this, we risk imprecision and
    // therefore glitches
    norm_type max_delta_norm = 0;

    // Returns the epsilon, and whether the epsilon is "accurate"
    // If false, then the first item is not set.
    std::pair<epsilon_type, bool> epsilon(delta_type delta) const {
      auto nd = fractals::norm(delta); // TODO: Avoid recomputing this
      if (nd > max_delta_norm)
        return {epsilon_type(), false};
      auto d = delta;
      epsilon_type s(0);
      typename HighExponentComplex::value_type prev_norm = 0, term_norm = 0;
      for (int t = 0; t < Terms; t++) {
        auto term = terms[t] * d;
        prev_norm = term_norm;
        term_norm = fractals::norm(term);
        s += term;
        d = d * delta;
      }
      return std::make_pair(s, true);
    }
  };

  taylor_series_orbit<HighExponentComplex, ReferenceOrbit, Terms> reference;
  std::vector<Entry> entries;

  void get_next() {
    entries.push_back({convert_complex<Complex>(*reference), reference.terms});
    ++reference;
  }

public:
  using relative_orbit =
      perturbation_orbit<value_type, value_type, stored_taylor_series_orbit>;

  // Returns a new relative orbit, with a certain number of iterations already
  // skipped.
  // - `limit` is the bailout value
  // - `iterations_skipped` is an in/out value for the number of iterations
  // already taken. It's used to optimize the search for the next iteration
  // limit.
  relative_orbit make_relative_orbit(delta_type delta, int limit,
                                     int &iterations_skipped) const {
    auto s = find_iterations_to_skip(delta, entries.size(), iterations_skipped);
    return {*this, convert_complex<Complex>(delta), iterations_skipped, s};
  }
};

} // namespace mandelbrot
