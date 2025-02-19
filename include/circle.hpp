/*
    A sample fractal showing how to add new fractals.

    The general idea is that you create a new symbol of type `const
    fractals::fractal &`. The implementation details of the
    fractal are hidden in the implementation file (circle.cpp)
*/

#include "mandelbrot_fwd.hpp"

extern const fractals::fractal &circle_fractal;
