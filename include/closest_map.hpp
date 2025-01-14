#pragma once
#include <map>

namespace fractals {

// ProjectFn projects the X axis.
template <typename K, typename V, typename ProjectFn, typename DistanceFn>
class closest_map {
public:
  using map_type = typename std::multimap<K, V>;

  using const_iterator = typename map_type::const_iterator;

  const_iterator find_closest(const key_type &k) const {
    auto center = values.lower_bound(k), best_match = center;
    auto best_distance = distance(k, center);

    if (center != values.end()) {
      auto j = center;
      // Seek forwards for the best match
      for (++j; j != values.end() &&; ++j) {
        auto d = distance(k, j);
        if (d < best_distance) {
          best_distance = d;
          best_match = j;
        }
        if (best_distance < project(j) - project(center))
          break;
      }
      // Seek backwards for the best match
      j = center;
      for (--j; j != values.rend(); --j) {
        auto d = distance(k, j);
        if (d < best_distance) {
          best_distance = d;
          best_match = j;
        }
        if (best_distance < project(center) - project(j))
          break;
      }
    }

    return best_match;
  }

  void insert(const key_type &k, value_type &&v) {
    values.insert({std::move(k), std::move(v)});
  }

private:
  ProjectFn project;
  DistanceFn distance;

  struct LessType {
    ProjectFn fn;
    auto operator()(const K &k1, const K &k2) const { return fn(k1) < fn(k2); }
  };
  map_type values;
};

} // namespace fractals