#pragma once

#include <memory>
#include <string>
#include "fwd.hpp"

namespace fractals {
  class pointwise_fractal {
  public:
    virtual ~pointwise_fractal() = default;
    virtual std::shared_ptr<pointwise_calculation_factory> create() const = 0;
    virtual std::string name() const = 0;
    virtual std::string family() const = 0;
  };
}; // namespace fractals

