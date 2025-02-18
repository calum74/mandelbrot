/*
    An orbit is an infinite sequence of complex numbers modelled as a C++
    forward iterator or random-access iterator.

    The first element of an orbit is 0.

    The current value of the orbit is given by operator*()
    The orbit is advanced using operator++().
    For random-access orbits, the value of an element is given by operator[].
*/

#pragma once
#include "mandelbrot_calculation.hpp"
#include <atomic>
#include <cassert>
#include <optional>
#include <vector>

namespace mandelbrot {

template <typename T>
concept IteratedOrbit = requires(T v) {
  *v;
  { *v } -> std::convertible_to<typename T::value_type>;
  ++v;
  // v.reset();
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
class naive_relative_orbit {
public:
  using value_type = OrbitType;
  using delta_type = DeltaType;
  using epsilon_type = delta_type;
  using reference_orbit_type = ReferenceOrbit;
  using calculation = typename reference_orbit_type::calculation;

  naive_relative_orbit(reference_orbit_type ref, delta_type delta)
      : reference_orbit(ref), delta{delta}, epsilon{} {}

  naive_relative_orbit(const naive_relative_orbit &other)
      : reference_orbit(other.reference_orbit), delta(other.delta), epsilon{} {}

  void reset() {
    reference_orbit.reset();
    epsilon = {};
  }

  value_type operator*() const {
    return *reference_orbit + convert<value_type>(epsilon);
  }

  value_type reference_z() const { return *reference_orbit; }

  naive_relative_orbit &operator++() {
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

  value_type operator[](int n) const {
    assert(n >= 0 && n < values.size());
    return values.at(n);
  }

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
naive_relative_orbit<typename Rel::value_type, Rel>
make_naive_relative_orbit(Rel rel, auto delta) {
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
  using term_type = typename fractals::normalized<TermType>::type;
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
    ++orbit;
    return *this;
  }

private:
  ReferenceOrbit orbit;
};

// An orbit that's relative to another reference orbit, so can be computed
// using a low-precision complex number. ReferenceOrbit must be a
// random-access orbit (supporting [])
template <Complex C, Complex DeltaType, RandomAccessOrbit ReferenceOrbit>
class perturbation_orbit {
public:
  using value_type = C;
  using calculation = typename ReferenceOrbit::calculation;
  using delta_type = DeltaType;
  using epsilon_type = delta_type;

  perturbation_orbit() = default; // Not great.

  perturbation_orbit(const ReferenceOrbit &ref, delta_type delta)
      : perturbation_orbit(ref, delta, 0, 0, {}) {}

  perturbation_orbit(const ReferenceOrbit &ref, delta_type delta,
                     int starting_iteration, int start_j,
                     epsilon_type starting_epsilon)
      : n{starting_iteration}, j{start_j}, delta(delta),
        epsilon{starting_epsilon}, reference{&ref} {
    if (j == ref.size() - 1) {
      epsilon = (*reference)[j] + convert<value_type>(epsilon);
      j = 0;
    }
  }

  // Constructs a perturbation orbit that has the same reference orbit
  // and same iteration number, but now refers to a different delta.
  perturbation_orbit
  split_relative(delta_type delta_from_current_orbit,
                 epsilon_type epsilon_from_current_orbit) const {
    return {*reference, delta + delta_from_current_orbit, n, j,
            epsilon_from_current_orbit + epsilon};
  }

  int iteration() const { return n; }

  value_type operator*() const {
    assert(j >= 0 && j < reference->size());
    return (*reference)[j] + convert<value_type>(epsilon);
  }

  perturbation_orbit &operator++() {

    epsilon = calculation::step_epsilon((*reference)[j], epsilon, delta);
    epsilon = fractals::normalize(epsilon);
    j++;

    assert(j >= 0 && j < reference->size());
    auto z = **this;
    if (j == (*reference).size() - 1 || escaped((*reference)[j]) ||
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

public:
  epsilon_type epsilon; // !! Bad encapsulation
  const ReferenceOrbit *reference = 0;
};

template <Complex C, Complex DeltaType, RandomAccessOrbit ReferenceOrbit>
perturbation_orbit<C, DeltaType, ReferenceOrbit>
make_perturbation_orbit(const ReferenceOrbit &ref, DeltaType delta) {
  return {ref, delta};
}

template <int Precision1 = 10, int Precision2 = 100, Complex TermType,
          unsigned long Terms>
typename TermType::value_type
maximum_delta_norm(const std::array<TermType, Terms> &terms) {
  using term_norm = typename TermType::value_type;

  term_norm max_delta_norm = 0;
  term_norm p1(Precision1);
  term_norm p2(Precision2);
  term_norm prev_norm = fractals::norm(terms[0]);

  for (int i = 1; i < Terms - 1; i++) {
    auto n = fractals::norm(terms[i]);
    auto nr = prev_norm / (p1 * n);
    prev_norm = n;
    if (i == 1 || max_delta_norm > nr)
      max_delta_norm = nr;
  }
  auto n = fractals::norm(terms[Terms - 1]);
  auto nr = prev_norm / (p2 * n);

  if (max_delta_norm > nr)
    max_delta_norm = nr;

  return max_delta_norm;
}

template <Complex DeltaType, Complex TermType, unsigned long Terms>
DeltaType evaluate_epsilon(DeltaType delta,
                           const std::array<TermType, Terms> &terms) {
  TermType d = convert<TermType>(delta);
  TermType dc = d;
  TermType result = terms[0] * dc;
  for (int i = 1; i < Terms; ++i) {
    d = d * dc;
    result += terms[i] * d;
  }
  return convert<DeltaType>(result);
}

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
          IteratedOrbit ReferenceOrbit, int Terms, int TermPrecision1,
          int TermPrecision2>
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

  auto epsilon(int i, delta_type delta, auto nd) const {
    assert(i >= 0 && i < entries.size());
    return entries.at(i).epsilon(delta, nd);
  }

  value_type operator[](int i) const {
    assert(i >= 0 && i < entries.size());
    return entries.at(i).z;
  }

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
    auto nd = fractals::norm(delta);

    int skipped = iterations_skipped;

    if (skipped > max)
      skipped = max;
    if (skipped >= this->entries.size())
      skipped = this->entries.size() - 1;

    auto e = this->epsilon(skipped, delta, nd);
    if (e && !escaped((*this)[skipped])) {
      // Seek upwards
      min = skipped;
      epsilon = convert<epsilon_type>(*e);
      for (int mid = skipped + window_size; mid < max;
           mid += (window_size *= 2)) {
        e = this->epsilon(mid, delta, nd);
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
      for (int mid = skipped - window_size; mid >= 0;
           mid -= (window_size *= 2)) {
        e = this->epsilon(mid, delta, nd);
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
      auto e = this->epsilon(mid, delta, nd);
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
    std::array<typename normalized<term_type>::type, Terms> terms;
    // The maximum delta for this term
    // If we use a delta with a higher norm than this, we risk imprecision and
    // therefore glitches
    delta_norm max_delta_norm = 0;

    Entry(value_type z, const std::array<typename normalized<term_type>::type, Terms> &ts)
        : z(z), terms(ts) {
      // We'll look at the terms in the series to figure out what the maximum
      // size of delta is for this term before imprecision sets in.
      // Each term has a "norm" giving an indication of its size.
      // We need to make sure that each term is sufficiently "small"
      // relative to the previous term.
      max_delta_norm = convert<delta_norm>(
          maximum_delta_norm<TermPrecision1, TermPrecision2>(terms));
    }

    // Returns the epsilon, if it's accurate.
    std::optional<epsilon_type> epsilon(delta_type delta, delta_norm nd) const {
      if (nd > max_delta_norm)
        return std::nullopt;
      return convert<epsilon_type>(
          evaluate_epsilon(convert<term_type>(delta), terms));
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
            iterations_skipped, convert<epsilon_type>(s)};
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
`std::complex<high_exponent_real>` depending on the precision required.

`ReferenceOrbit` is a high precision reference orbit, that is iterated
precisely once.
*/
template <Complex OrbitType, Complex DeltaType, Complex TermType,
          IteratedOrbit ReferenceOrbit, int Terms, int TermPrecision1,
          int TermPrecision2>
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
      perturbation_orbit<value_type, delta_type, primary_orbit_type>;

  using secondary_orbit_type =
      stored_taylor_series_orbit<value_type, delta_type, term_type,
                                 secondary_reference_type, Terms,
                                 TermPrecision1, TermPrecision2>;

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
    return {{primary_orbit, delta}, max_iterations, stop};
  }
};

} // namespace mandelbrot
