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

template <typename T>
concept IteratedOrbit = requires(T v) {
  *v;
  { *v } -> std::convertible_to<typename T::value_type>;
  ++v;
  v.reset();
  // requires T::value_type;
  // requires(T::value_type);
  // requires(T::calculation);
};

template <typename T>
concept RandomAccessOrbit = requires(T v, int i) {
  { v[i] } -> std::same_as<typename T::value_type>;
};

template <typename T>
concept Calculation =
    requires(T) { T::step(std::complex<double>(), std::complex<double>()); };

/*
  A basic_orbit iterates the escape sequence manually.
  Often this is used as a "reference orbit" for more advanced algorithms.

  Complex - the type of the complex number
  Calculation - how to calculate the sequence
*/
template <Complex T, Calculation C> class basic_orbit {
public:
  using value_type = T;
  using calculation = C;

  basic_orbit(value_type c = value_type{}) : c{c}, z{} {}

  void reset() { z = 0; }

  const value_type &operator*() const { return z; }

  basic_orbit &operator++() {
    z = calculation::step(z, c);
    return *this;
  }

private:
  value_type c, z;
};

template <Calculation C, Complex T>
basic_orbit<T, C> make_basic_orbit(const T &c) {
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
template <Complex T, IteratedOrbit ReferenceOrbit> class converted_orbit {
public:
  converted_orbit() = default;
  converted_orbit(const ReferenceOrbit &r) : reference(r) {}

  using value_type = T;
  using calculation = typename ReferenceOrbit::calculation;

  value_type operator*() const { return convert<value_type>(*reference); }

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
template <Complex LowPrecisionComplex, Complex HighPrecisionComplex,
          Calculation C>
using high_precision_orbit =
    converted_orbit<LowPrecisionComplex, basic_orbit<HighPrecisionComplex, C>>;

template <Complex OrbitType, IteratedOrbit ReferenceOrbit,
          Complex DeltaType = OrbitType>
class relative_orbit {
public:
  using value_type = OrbitType;
  using delta_type = DeltaType;
  using epsilon_type = delta_type;
  using reference_orbit_type = ReferenceOrbit;
  using calculation = typename reference_orbit_type::calculation;

  relative_orbit(reference_orbit_type ref, value_type delta)
      : reference_orbit(ref), delta{delta}, epsilon{} {}

  relative_orbit(const relative_orbit &other)
      : reference_orbit(other.reference_orbit), delta(other.delta), epsilon{} {}

  void reset() {
    reference_orbit.reset();
    epsilon = {};
  }

  value_type operator*() const { return *reference_orbit + epsilon; }

  value_type reference_z() const { return *reference_orbit; }

  relative_orbit &operator++() {
    epsilon = calculation::step_epsilon(*reference_orbit, epsilon, delta);
    ++reference_orbit;
    return *this;
  }

  value_type get_epsilon() const { return epsilon; }

private:
  reference_orbit_type reference_orbit;
  delta_type delta;
  epsilon_type epsilon;
};

template <Complex C, IteratedOrbit Ref> class stored_orbit {
public:
  stored_orbit(Ref orbit, int max_iterations, std::atomic<bool> &stop) {
    // Store up to max_iterations or until orbit has escaped
    // Note that `orbit` is a value type so is reset to the beginning
    do {
      values.push_back(convert<C>(*orbit));
      ++orbit;
    } while (!stop && values.size() <= max_iterations &&
             !escaped(values.back()));
  }

  using value_type = C;
  using calculation = typename Ref::calculation;

  value_type operator[](int n) const { return values.at(n); }

  struct reference_orbit {
  public:
    using value_type = C;
    using calculation = typename Ref::calculation;
    reference_orbit(const stored_orbit &ref) : ref{ref}, current{0} {}
    reference_orbit(const reference_orbit &ref) : ref{ref.ref}, current{0} {}

    value_type operator*() const { return ref[current]; }

    reference_orbit &operator++() {
      ++current;
      return *this;
    }

    void reset() { current = 0; }

  private:
    int current;

    const stored_orbit &ref;
  };

  reference_orbit make_reference() const { return {*this}; }

  auto size() const { return values.size(); }

private:
  std::vector<value_type> values;
};

template <Complex C, IteratedOrbit Ref>
stored_orbit<C, Ref> make_stored_orbit(const Ref &r, int max_iterations,
                                       std::atomic<bool> &stop) {
  return {r, max_iterations, stop};
}

template <IteratedOrbit Rel>
relative_orbit<typename Rel::value_type, Rel> make_relative_orbit(Rel rel,
                                                                  auto delta) {
  return {rel, delta};
}

/*
  An orbit that also computes terms of a Taylor series approximation for
  each iteration.

  z = z_0 + A.delta + B.delta^2 + C.delta^3
*/
template <Complex OrbitType, Complex TermType, IteratedOrbit ReferenceOrbit,
          int Terms>
class taylor_series_orbit {
public:
  using value_type = OrbitType;
  using term_type = TermType;
  using calculation = typename ReferenceOrbit::calculation;

  std::array<term_type, Terms> terms;

  taylor_series_orbit() = default;

  taylor_series_orbit(ReferenceOrbit o) : orbit(o) {}

  // Convert the high-precision value into a low-precision value
  value_type operator*() const { return convert<value_type>(*orbit); }

  taylor_series_orbit &operator++() {

    // See
    // https://en.wikipedia.org/wiki/Plotting_algorithms_for_the_Mandelbrot_set
    terms = calculation::delta_terms(**this, terms);
    for (auto &t : terms)
      t = normalize(t);
    ++orbit;
    return *this;
  }

private:
  ReferenceOrbit orbit;
};

// An orbit that's relative to another reference orbit, so can be computed
// using a low-precision complex number. ReferenceOrbit must be a
// random-access orbit (supporting [])
template <Complex C, Complex HighExponentComplex,
          RandomAccessOrbit ReferenceOrbit>
class perturbation_orbit {
public:
  using value_type = C;
  using calculation = typename ReferenceOrbit::calculation;
  using delta_type = HighExponentComplex;
  using epsilon_type = delta_type;

  perturbation_orbit(const ReferenceOrbit &ref, delta_type delta,
                     int starting_iteration = 0,
                     epsilon_type starting_epsilon = {})
      : n{starting_iteration}, j{n}, delta(delta), epsilon{starting_epsilon},
        reference{ref} {}

  perturbation_orbit(const perturbation_orbit &other) = delete;

  int iteration() const { return n; }

  value_type operator*() const {
    return reference[j] + convert<value_type>(epsilon);
  }

  perturbation_orbit &operator++() {

    epsilon = calculation::step_epsilon(reference[j], epsilon, delta);
    epsilon = fractals::normalize(epsilon);
    j++;

    auto z = reference[j] + convert<value_type>(epsilon);

    if (escaped(reference[j]) ||
        fractals::norm(z) < fractals::norm(convert<value_type>(epsilon))) {
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
template <Complex OrbitType, Complex DeltaType, Complex TermType,
          IteratedOrbit ReferenceOrbit, int Terms, int Precision>
class stored_taylor_series_orbit {
public:
  using value_type = OrbitType;
  using term_type = TermType;
  using delta_type = DeltaType;
  using epsilon_type = DeltaType;
  using calculation = typename ReferenceOrbit::calculation;

  stored_taylor_series_orbit() = default;
  stored_taylor_series_orbit(ReferenceOrbit r) : reference(r) {}
  stored_taylor_series_orbit(ReferenceOrbit r, int max_iterations,
                             std::atomic<bool> &stop)
      : reference(r) {
    do {
      get_next();
    } while (!stop && entries.size() <= max_iterations &&
             !escaped(entries.back().z));
    // debug_terms();
  }

  // !! This does not belong here
  void debug_terms() const {
    // Debug - look at the terms
    int max_term[Terms] = {};
    using R = TermType::value_type;
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

  auto epsilon(int i, delta_type delta) const {
    return entries.at(i).epsilon(delta);
  }

  value_type operator[](int i) const { return entries.at(i).z; }

  auto size() const { return entries.size(); }

private:
  // Finds the number of iterations it's safe to skip
  // because we judge that we are sufficiently close to the reference orbit
  epsilon_type find_iterations_to_skip(delta_type delta, int max,
                                       int &iterations_skipped) const {
    int min = 0;

    // Step 1: Establish the max and min
    epsilon_type epsilon = {};
    int window_size = 4;

    int skipped = iterations_skipped;

    if (skipped > max)
      skipped = max;
    if (skipped >= this->entries.size())
      skipped = this->entries.size() - 1;

    auto e = this->epsilon(skipped, delta);
    if (e && !escaped((*this)[skipped])) {
      // Seek upwards
      min = skipped;
      epsilon = convert<epsilon_type>(*e);
      for (int mid = iterations_skipped + window_size; mid < max;
           mid += (window_size *= 2)) {
        e = this->epsilon(mid, delta);
        if (e && !escaped((*this)[mid])) {
          min = mid;
          epsilon = *e;
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
        if (e && !escaped((*this)[mid])) {
          min = mid;
          epsilon = *e;
          break;
        } else {
          max = mid;
        }
      }
    }

    // Step 2: Find the min using binary search
    while (max - min > 4) {
      int mid = (max + min) / 2;
      auto e = this->epsilon(mid, delta);
      if (e && !escaped((*this)[mid])) {
        min = mid;
        epsilon = *e;
      } else
        max = mid;
    }
    iterations_skipped = min;

    return epsilon;
  }

private:
  // An `Entry` is an orbital value (`z`) and the Taylor series terms
  // (`terms`) for an iteration.
  struct Entry {

    using delta_norm = typename delta_type::value_type;
    using term_norm = typename term_type::value_type;

    value_type z;
    std::array<term_type, Terms> terms;
    // The maximum delta for this term
    // If we use a delta with a higher norm than this, we risk imprecision and
    // therefore glitches
    delta_norm max_delta_norm = 0;

    Entry(value_type z, const std::array<term_type, Terms> &ts)
        : z(z), terms(ts) {
      // We'll look at the terms in the series to figure out what the maximum
      // size of delta is for this term before imprecision sets in.
      // Each term has a "norm" giving an indication of its size.
      // We need to make sure that each term is sufficiently "small"
      // relative to the previous term.

      auto prev_norm = fractals::norm(terms[0]);
      for (int i = 1; i < Terms - 1; i++) {
        auto n = fractals::norm(terms[i]);
        auto nr = convert<delta_norm>(prev_norm /
                                      (term_norm(10) * n)); // TODO: Tweak this
        prev_norm = n;
        if (i == 1 || max_delta_norm > nr)
          max_delta_norm = nr;
      }
      auto n = fractals::norm(terms[Terms - 1]);
      auto nr = convert<delta_norm>(
          prev_norm /
          (term_norm(100) * n)); // TODO: Tweak this - Precision is ignored

      if (max_delta_norm > convert<delta_norm>(nr))
        max_delta_norm = nr;
    }

    // Returns the epsilon, if it's accurate.
    std::optional<epsilon_type> epsilon(delta_type delta) const {
      auto nd = fractals::norm(delta); // TODO: Avoid recomputing this
      if (nd > max_delta_norm)
        return std::nullopt;
      auto d_conv = convert<term_type>(delta);
      term_type d = d_conv;
      term_type s(0);
      typename TermType::value_type prev_norm = 0, term_norm = 0;
      for (int t = 0; t < Terms; t++) {
        auto term = terms[t] * d;
        prev_norm = term_norm;
        term_norm = fractals::norm(term);
        s += term;
        d = d * d_conv;
      }
      return convert<epsilon_type>(s);
    }
  };

  taylor_series_orbit<value_type, term_type, ReferenceOrbit, Terms> reference;
  std::vector<Entry> entries;

  void get_next() {
    entries.push_back({convert<value_type>(*reference), reference.terms});
    ++reference;
  }

public:
  using relative_orbit =
      perturbation_orbit<value_type, epsilon_type, stored_taylor_series_orbit>;

  // Returns a new relative orbit, with a certain number of iterations already
  // skipped.
  // - `limit` is the bailout value
  // - `iterations_skipped` is an in/out value for the number of iterations
  // already taken. It's used to optimize the search for the next iteration
  // limit.
  relative_orbit make_relative_orbit(delta_type delta, int limit,
                                     int &iterations_skipped) const {
    auto s = find_iterations_to_skip(delta, entries.size(), iterations_skipped);
    return {*this, convert<epsilon_type>(delta), iterations_skipped,
            convert<epsilon_type>(s)};
  }
};

/*
A group of "secondary" reference orbits clustered around a single
high precision "primary" reference orbit. The big idea is to make it cheap to
create secondary orbits that are closer to the point being calculated,
allowing it to skip a much higher number of iterations.

`OrbitType` is a regular complex number like `std::complex<double>`.

`TermType` is used to represent terms in the Taylor series. These should be
high exponent values such as `std::complex<high_exponent_real>`. Otherwise the
Taylor series will run out of iterations even for low precisions.

`DeltaType` is used to represent deltas, and can be `std::complex<double>` or
`std::complex<high_exponent_real>` depending on the current precision being
calculated.

`ReferenceOrbit` is a high precision reference orbit, that is iterated
precisely once.
*/
template <Complex OrbitType, Complex DeltaType, Complex TermType,
          IteratedOrbit ReferenceOrbit, int Terms, int TermPrecision>
class taylor_series_cluster {
public:
  using calculation = typename ReferenceOrbit::calculation;
  using value_type = OrbitType;
  using delta_type = DeltaType;
  using term_type = TermType;
  using epsilon_type = DeltaType;

  // We only need to store the primary orbit in low precision
  // (although obviously it was initially calculated to high precision)
  using primary_orbit_type = stored_orbit<value_type, ReferenceOrbit>;
  primary_orbit_type primary_orbit;

  using secondary_reference_type =
      relative_orbit<value_type, typename primary_orbit_type::reference_orbit,
                     delta_type>;

  using secondary_orbit_type =
      stored_taylor_series_orbit<value_type, delta_type, term_type,
                                 secondary_reference_type, Terms,
                                 TermPrecision>;

  using relative_orbit = typename secondary_orbit_type::relative_orbit;

  // We initially calculate the primary orbit, but don't create any secondary
  // orbits yet
  taylor_series_cluster(ReferenceOrbit orbit, int max_iterations,
                        std::atomic<bool> &stop)
      : primary_orbit(orbit, max_iterations, stop) {}

  // Constructs and fully evaluates a secondary orbit
  // (thread-safe)
  secondary_orbit_type make_secondary_orbit(delta_type delta,
                                            int max_iterations,
                                            std::atomic<bool> &stop) const {
    return {{primary_orbit.make_reference(), delta}, max_iterations, stop};
  }

  secondary_orbit_type get_closest_orbit(delta_type);

  relative_orbit make_relative_orbit(delta_type delta, int limit,
                                     int &iterations_skipped) const;

  //
  // secondary_orbit_type central_orbit;

  // Our list of
  std::vector<secondary_orbit_type> secondary_orbits;
};
} // namespace mandelbrot
