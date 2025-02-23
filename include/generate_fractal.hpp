#pragma once

#include "make_fractal.hpp"

namespace fractals {
    template <typename Fractal>
    detail::MultiPrecisionFractal<
        typename Fractal::template precision<50>::type, typename Fractal::template precision<128>::type,
        typename Fractal::template precision<256>::type, typename Fractal::template precision<640>::type,
        typename Fractal::template precision<1024>::type, typename Fractal::template precision<1596>::type,
        typename Fractal::template precision<2048>::type, typename Fractal::template precision<4096>::type>
    generate_fractal(const char *name, const char * family) {
      return {name, family};
    }
    
    template <typename Fractal>
    auto generate_fractal(const char *name) {
      return generate_fractal<Fractal>(name, name);
    }    
}