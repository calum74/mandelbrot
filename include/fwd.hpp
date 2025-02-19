#pragma once

namespace fractals
{
class pointwise_calculation_factory;
class pointwise_fractal;
class pointwise_calculation;
class view_coords;

template <typename T1, typename T2> T1 convert(const T2 &x);
template <typename C1, typename C2=C1> class plane;
}