/*
    A sample fractal showing how to add new fractals.

    The general idea is that you create a new symbol of type `const
    fractals::pointwise_fractal &`. The implementation details of the
    fractal are hidden in the implementation file (circle.cpp)
*/

#include "fwd.hpp"

extern const fractals::pointwise_fractal &circle_fractal;
