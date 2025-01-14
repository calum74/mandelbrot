#pragma once
#include <map>

namespace fractals {

// ProjectFn projects the X axis.
template <typename K, typename V, typename ProjectFn, typename DistanceFn>
class closest_map {
public:
  struct LessType {
    ProjectFn fn;
    bool operator()(const K &k1, const K &k2) const { return fn(k1) < fn(k2); }
  };

  using map_type = typename std::multimap<K, V, LessType>;

  using const_iterator = typename map_type::const_iterator;

  const_iterator find_closest(const K &k) const {
    if (values.begin() == values.end())
      return values.end();

    auto center = values.lower_bound(k);
    if (center == values.end())
      center--;
    auto best_match = center;

    auto best_distance = distance(k, center->first);
    auto j = center;
    // Seek forwards for the best match
    for (++j; j != values.end(); ++j) {
      auto d = distance(k, j->first);
      if (d < best_distance) {
        best_distance = d;
        best_match = j;
      }
      if (best_distance < project(j->first) - project(center->first))
        break;
    }

    if (center != values.begin()) {
      // Seek backwards for the best match
      auto j = center;
      for (--j;; --j) {
        auto d = distance(k, j->first);
        if (d < best_distance) {
          best_distance = d;
          best_match = j;
        }
        if (best_distance < project(center->first) - project(j->first))
          break;
        if (j == values.begin())
          break;
      }
    }

    return best_match;
  }

  void insert(const K &k, V &&v) {
    values.insert({std::move(k), std::move(v)});
  }

private:
  ProjectFn project;
  DistanceFn distance;
  map_type values;
};

} // namespace fractals