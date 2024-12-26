/*
    A sample fractal showing how to add new fractals.

    The general idea is that you create a new symbol of type `const
    fractals::PointwiseFractal &`. The implementation details of the fractal are
    hidden in the implementation file (circle.cpp)
*/

#include "fractal.hpp"

extern const fractals::PointwiseFractal &circle_fractal;
