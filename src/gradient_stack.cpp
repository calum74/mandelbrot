#include "gradient_stack.hpp"

fractals::gradient_stack::result
fractals::gradient_stack::map_iteration(double d, double default_gradient,
                                        double default_offset) const {
  for (auto j = stack.rbegin(); j != stack.rend(); ++j) {
    if (d > j->iteration) {
      return {d / j->gradient + j->offset, j->gradient};
    }
  }

  return {d / default_gradient + default_offset, default_gradient};
}

void fractals::gradient_stack::push(double iteration, double new_gradient,
                                    double default_gradient,
                                    double default_offset) {
  // Remove any colours that are above the current max
  while (!stack.empty() && stack.back().iteration > iteration) {
    stack.pop_back();
  }

  auto last_gradient = stack.empty() ? default_gradient : stack.back().gradient;
  auto last_offset = stack.empty() ? default_offset : stack.back().offset;

  /*
  To align the colours, we need to ensure that
  iteration/last_gradient + last_offset = iteration/new_gradient + new_offset
  -> new_offset = last_offset + iteration/last_gradient -
  iteration/new_gradient
  */

  auto new_offset =
      last_offset + iteration / last_gradient - iteration / new_gradient;
  stack.push_back({iteration, new_gradient, new_offset});
}

void fractals::gradient_stack::clear() { stack.clear(); }