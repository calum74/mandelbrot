cmake_minimum_required(VERSION 3.21)
project(mandelbrot)
set(CMAKE_CXX_STANDARD 23)
enable_testing()

add_library(mandelbrot complex.hpp convert.hpp fractal.hpp high_precision_real.hpp mandelbrot_calculation.hpp 
    mandelbrot.hpp plane.hpp orbit.hpp layer.cpp layer.hpp mandelbrot.cpp mandeldrop.cpp view_coords.hpp view_coords.cpp rendering_sequence.cpp rendering_sequence.hpp)
add_library(extra_fractals circle.cpp rendering_sequence.cpp rendering_sequence.hpp)
link_libraries(mandelbrot)

add_executable(tests tests.cpp)
add_executable(high_precision_tests high_precision_tests.cpp)
add_executable(tutorial tutorial.cpp)

add_test(unittests tests)
add_test(hp high_precision_tests)
add_test(tutorial tutorial)
