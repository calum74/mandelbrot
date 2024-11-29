#pragma once

#include "fractal.hpp"

extern const fractals::PointwiseFractal &mb;

namespace mandelbrot {
void add_fractals(fractals::Registry &);
}
