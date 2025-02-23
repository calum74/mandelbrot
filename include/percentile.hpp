#pragma once
#include <algorithm>

namespace fractals {
// Helper to find high percentiles (e.g. 99.9%)
// Technically O(n.ln n) but faster when we're looking for high percentiles
// because only the top items are sorted.
template <typename It, typename Cmp> It top_percentile(It first, It last, double p, Cmp cmp = std::less<typename It::value_type>()) {
  std::make_heap(first, last, cmp);
  int to_skip = std::distance(first, last) * (1.0 - p);
  for (int i = 0; i <= to_skip; ++i, --last)
    std::pop_heap(first, last, cmp);
  return last;
}
} // namespace util