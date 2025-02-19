#pragma once

#include <memory>

#include "fwd.hpp"
#include "plane.hpp"
#include "view_coords.hpp"

namespace fractals {
  
  class pointwise_fractal {
  public:
    virtual ~pointwise_fractal() = default;
    virtual std::shared_ptr<pointwise_calculation_factory> create() const = 0;
    virtual std::string name() const = 0;
    virtual std::string family() const = 0;
  };
  



}; // namespace fractals

