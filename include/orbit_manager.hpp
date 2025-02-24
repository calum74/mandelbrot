#pragma once
#include "orbit.hpp"
#include "rendering_sequence.hpp"
#include <future>
#include <memory>

namespace mandelbrot {

/*
The orbit manager is used to compute a group of reference orbits.

In order to keep the UI responsive, it computes reference orbits in the
background, and automatically selects the best reference orbit. As more orbits
are computed, so better and better reference orbits are returned.

A "better" reference orbit is closer to the point being computed, so is likely
to be able to skip more iterations.
*/
template <Complex OrbitType, Complex DeltaType, Complex TermType,
          unsigned long Terms, int TermPrecision1, int TermPrecision2,
          IteratedOrbit HighPrecisionReferenceOrbit>
class orbit_manager {

  using cluster_type = taylor_series_cluster<OrbitType, DeltaType, TermType,
                                             HighPrecisionReferenceOrbit, Terms,
                                             TermPrecision1, TermPrecision2>;

  using primary_orbit_type = typename cluster_type::primary_orbit_type;

  using secondary_orbit_type = typename cluster_type::secondary_orbit_type;
  using relative_orbit = typename cluster_type::relative_orbit;

  struct secondary_orbit {
    secondary_orbit(DeltaType delta, secondary_orbit_type &&orbit)
        : delta(delta), orbit(std::move(orbit)) {}
    DeltaType delta;
    secondary_orbit_type orbit;
    std::atomic<int> iterations_skipped_cache = 0;
  };

  using secondary_reference_type =
      typename cluster_type::secondary_reference_type;

private:
  // Seeds the orbit_manager completely from scratch.
  // In this case, the entire reference orbit needs to be computed
  // and the series for this orbit. This could be a bit slow.
  // Not threadsafe
  void initialize(const HighPrecisionReferenceOrbit &init, int max_iterations,
                  std::atomic<bool> &stop) {

    primary_orbit_type new_primary(init, max_iterations, stop);

    if (stop)
      return;

    auto series = std::make_shared<secondary_orbit>(
        DeltaType{},
        secondary_orbit_type{secondary_reference_type{new_primary, DeltaType{}},
                             max_iterations, stop});
    if (stop)
      return;

    orbit_storage.clear();
    primary_series = std::move(series);
  }

public:
  // Starts the orbit manager by reusing the current reference orbit. This
  // is fast. delta is the delta from the previous center. Not threadsafe
  void new_view(DeltaType delta, DeltaType maxDelta, int maxSecondaryOrbits,
                const HighPrecisionReferenceOrbit &init, int max_iterations,
                std::atomic<bool> &stop) {

    try_stop_orbits_thread();

    // Determine whether the previous reference orbit is suitable as a temporary
    // reference orbit. Heuristically, the orbit
    // - must be within the current view
    // - has escaped, or
    // - has not escaped, and is sufficiently long
    if (!primary_series ||
        fractals::norm(primary_series->delta - delta) >
            fractals::norm(maxDelta) ||
        (!escaped(primary_series->orbit[primary_series->orbit.size() - 1]) &&
         primary_series->orbit.size() < max_iterations)) {
      // Unfortunately, recycling the primary orbit isn't always feasible if
      // std::cout << "Info: synchronously recalculating primary orbit\n";
      // Actually calculate the primary orbit extra long to increase probability
      // of reuse
      initialize(init, 2 * max_iterations, stop);
      if (stop)
        return;
      delta = DeltaType{0};
    }

    orbit_storage.clear();
    orbit_storage.push_back(primary_series);

    auto new_primary_series = std::make_shared<secondary_orbit>(
        primary_series->delta - delta, std::move(primary_series->orbit));
    orbit_storage.push_back(new_primary_series);

    std::vector<std::atomic<secondary_orbit *>> new_lookup(
        maxSecondaryOrbits == 0 ? 1 : maxSecondaryOrbits * maxSecondaryOrbits);

    primary_series = new_primary_series;
    std::fill(new_lookup.begin(), new_lookup.end(), new_primary_series.get());

    orbit_lookup.swap(new_lookup);

    this->max_delta = maxDelta;
    lookup_width = maxSecondaryOrbits;
    lookup_height = maxSecondaryOrbits;

    orbits_thread = std::async([init, max_iterations, this] {
      thread_fn(init, max_iterations, stop_orbits_thread);
    });
  }

  ~orbit_manager() { try_stop_orbits_thread(); }

private:
  void try_stop_orbits_thread() {
    if (orbits_thread.valid()) {
      stop_orbits_thread = true;
      orbits_thread.wait();
      stop_orbits_thread = false;
    }
  }

  std::atomic<bool> stop_orbits_thread;
  std::future<void> orbits_thread;

  // This is a "slow" function that is intended to be in a thread and
  // gradually populates the orbits. In the meantime, the orbit manager is
  // fully usable but may return a slightly less optimals orbits (=slower)
  // until the orbits are fully populated.
  //
  // init is the true reference orbit at the center of the image
  // maxDelta is the maximum delta we are expecting (in the positive
  // direction) maxSecondary orbits maximum number of secondary orbits we
  // should compute.
  //
  // For now, only one thread function is available.
  // Threadsafe in one thread
  void thread_fn(HighPrecisionReferenceOrbit init, int max_iterations,
                 std::atomic<bool> &stop) {

    // 1) Compute new primary reference orbit
    // Make it extra long (*2) so we can recycle it in the next view
    primary_orbit_type new_primary(init, max_iterations * 2, stop);

    if (stop)
      return;

    auto series = std::make_shared<secondary_orbit>(
        DeltaType{},
        secondary_orbit_type{secondary_reference_type{new_primary, DeltaType{}},
                             max_iterations, stop});
    if (stop)
      return;

    primary_series = series;
    orbit_storage.push_back(series);
    std::fill(orbit_lookup.begin(), orbit_lookup.end(), primary_series.get());

    // 2) Populate secondary reference orbits

    if (lookup_width == 0 || lookup_height == 0)
      return;

    rendering_sequence rs(lookup_width, lookup_height, 16);
    int x, y, stride;
    bool stride_changed;
    while (!stop && rs.next(x, y, stride, stride_changed)) {
      // Construct a new orbit
      using DeltaReal = typename DeltaType::value_type;
      DeltaType orbit_delta{
          max_delta.real() *
              number_cast<DeltaReal>(double(x + 0.5) / double(lookup_width) - 0.5),
          max_delta.imag() *
              number_cast<DeltaReal>(double(y + 0.5) / double(lookup_height) -
                                 0.5)};

      auto series = std::make_shared<secondary_orbit>(
          orbit_delta, secondary_orbit_type{
                           secondary_reference_type{new_primary, orbit_delta},
                           max_iterations, stop});
      orbit_storage.push_back(series);

      for (int j = 0; j < stride && y + j < lookup_height; ++j)
        for (int i = 0; i < stride && i + x < lookup_width; ++i) {
          orbit_lookup[x + i + (y + j) * lookup_width] = series.get();
        }
    }
  }

public:
  // Returns the orbit for the point delta (relative to the center)
  // The orbit has already been fast-forwarded.
  // Threadsafe
  relative_orbit lookup(DeltaType delta, int max_iterations) const {
    // Find the cell corresponding to delta
    int x = 0.5 * lookup_width *
            (1.0 + number_cast<double>(delta.real() / max_delta.real()));
    int y = 0.5 * lookup_height *
            (1.0 + number_cast<double>(delta.imag() / max_delta.imag()));

    if (x >= lookup_width)
      x = lookup_width - 1;
    if (x < 0)
      x = 0;
    if (y >= lookup_height)
      y = lookup_height - 1;
    if (y < 0)
      y = 0;

    // Threadsafe because we are doing an atomic read of the orbit_lookup.
    // The underlying objects are guaranteed to be valid.
    auto index = x + y * lookup_width;
    assert(index >= 0);
    assert(index < orbit_lookup.size());
    secondary_orbit *local_reference = orbit_lookup.at(x + y * lookup_width);

    int skipped = local_reference->iterations_skipped_cache;
    auto result = local_reference->orbit.make_relative_orbit(
        delta - local_reference->delta, max_iterations, skipped);
    local_reference->iterations_skipped_cache = skipped;
    return result;
  }

private:
  std::shared_ptr<secondary_orbit> primary_series;

  // How orbits are indexed
  // Orbits are stored in a grid.
  // Deltas are mapped to grid positions
  // As secondary orbits are added, they replace orbits in the grid
  std::vector<std::shared_ptr<secondary_orbit>> orbit_storage;
  std::vector<std::atomic<secondary_orbit *>> orbit_lookup;

  int lookup_width, lookup_height;
  DeltaType max_delta;
};

} // namespace mandelbrot