#pragma once

#include <memory>
#include <string>
#include "mandelbrot_fwd.hpp"

namespace fractals {
  class fractal {
  public:
    virtual ~fractal() = default;
    virtual std::shared_ptr<fractal_calculation_factory> create() const = 0;
    virtual std::string name() const = 0;
    virtual std::string family() const = 0;
  };
}; // namespace fractals

