#include "fractal.h"
#include "mandelbrot.hpp"
#undef NDEBUG
#include <cassert>

int main()
{
    using namespace mandelbrot;
    using C = complex<double>;

    assert((C{0,0} + C{1,1} == C{1,1}));
    assert((square(C{0,1}) == C{-1,0}));

    return 0;
}
