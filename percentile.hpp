#pragma once
#include <algorithm>

namespace util {
// Helper to find the 99.9th percentile
// Technically O(n.ln n) but faster because we only sort the top 0.1% of
// the data
template <typename It>
typename It::value_type percentile(It first, It last, double p) {
  std::make_heap(first, last);
  int to_skip = std::distance(first, last) * (1.0 - p);
  for (int i = 0; i < to_skip; ++i, --last)
    std::pop_heap(first, last);
  return *last;
}
} // namespace util