#pragma once
#include "closest_map.hpp"
#include "orbit.hpp"

namespace mandelbrot {

// Manages a set of secondary reference orbits.
template <Complex OrbitType, Complex DeltaType, Complex TermType,
          IteratedOrbit ReferenceOrbit, int Terms, int TermPrecision>
class reference_orbit_manager
    : public taylor_series_cluster<OrbitType, DeltaType, TermType,
                                   ReferenceOrbit, Terms, TermPrecision> {

  using base_type = taylor_series_cluster<OrbitType, DeltaType, TermType,
                                          ReferenceOrbit, Terms, TermPrecision>;

public:
  using secondary_orbit_type = base_type::secondary_orbit_type;

  // You can add a secondary reference orbit at
  // any time. The orbit was obtained from
  // make_secondary_orbit(). Not threadsafe.
  void add_secondary_reference_orbit(DeltaType delta,
                                     secondary_orbit_type &&secondary);

  using map_type = closest_map<DeltaType, secondary_orbit_type>;

  // We don't know which secondary reference
  // orbit will be the best, but we choose the
  // closest. This delta is relative to the
  // primary reference orbit. When we query the
  // secondary orbit, its delta is relative to
  // the secondary orbit, so you might need to
  // subtract the primary delta. Threadsafe for
  // reading.
  typename map_type::iterator
  locate_closest_secondary_reference_orbit(DeltaType delta) const;

private:
  map_type orbits;
};

} // namespace mandelbrot
