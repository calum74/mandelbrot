#pragma once
#include <vector>

namespace fractals {
class gradient_stack {
public:
  struct result {
    double value;
    double gradient;
  };

  result map_iteration(double iteration, double default_gradient,
                       double default_offset) const;
  void push(double iteration, double new_gradient, double default_gradient,
            double default_offset);
  void clear();

private:
  struct entry {
    double iteration;
    double gradient;
    double offset;
  };
  entry default_entry;
  std::vector<entry> stack;
};
} // namespace fractals