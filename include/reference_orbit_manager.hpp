#pragma once
#include "closest_map.hpp"
#include "orbit.hpp"

namespace mandelbrot {

namespace detail {
template <Complex DeltaType> struct map_project {
  auto operator()(const DeltaType &x) const { return x.real() * x.real(); }
};

template <Complex DeltaType> struct map_distance {
  auto operator()(const DeltaType &x, const DeltaType &y) const {
    return fractals::norm(x - y);
  }
};

} // namespace detail

// Manages a set of secondary reference orbits.
template <Complex OrbitType, Complex DeltaType, Complex TermType,
          IteratedOrbit ReferenceOrbit, int Terms, int TermPrecision>
class reference_orbit_manager
    : public taylor_series_cluster<OrbitType, DeltaType, TermType,
                                   ReferenceOrbit, Terms, TermPrecision> {

  using base_type = taylor_series_cluster<OrbitType, DeltaType, TermType,
                                          ReferenceOrbit, Terms, TermPrecision>;

public:
  reference_orbit_manager(ReferenceOrbit orbit, int max_iterations,
                          std::atomic<bool> &stop)
      : base_type(orbit, max_iterations, stop) {}

  using secondary_orbit_type = base_type::secondary_orbit_type;

  // You can add a secondary reference orbit at
  // any time. The orbit was obtained from
  // make_secondary_orbit(). Not threadsafe.
  void add_secondary_reference_orbit(DeltaType delta,
                                     secondary_orbit_type &&secondary) {
    orbits.insert(delta, std::move(secondary));
  }

  using map_type = closest_map<DeltaType, secondary_orbit_type,
                               detail::map_project<DeltaType>,
                               detail::map_distance<DeltaType>>;

  // We don't know which secondary reference
  // orbit will be the best, but we choose the
  // closest. This delta is relative to the
  // primary reference orbit. When we query the
  // secondary orbit, its delta is relative to
  // the secondary orbit, so you might need to
  // subtract the primary delta. Threadsafe for
  // reading.
  typename map_type::const_iterator
  locate_closest_secondary_reference_orbit(DeltaType delta) const {
    return orbits.find_closest(delta);
  }

private:
  map_type orbits;
};

} // namespace mandelbrot
