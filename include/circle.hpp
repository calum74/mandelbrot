/*
    A sample fractal showing how to add new fractals.

    The general idea is that you create a new symbol of type `const
    fractals::pointwise_calculationFactory &`. The implementation details of the
   fractal are hidden in the implementation file (circle.cpp)
*/

#include "fractal.hpp"

extern const fractals::pointwise_calculationFactory &circle_fractal;
