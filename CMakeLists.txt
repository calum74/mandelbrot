cmake_minimum_required(VERSION 3.21)
project(mandelbrot)
set(CMAKE_CXX_STANDARD 23)
enable_testing()

add_library(mandelbrot mandelbrot.cpp)
add_library(extra_fractals circle.cpp)
link_libraries(mandelbrot)

add_executable(tests tests.cpp)
add_executable(high_precision_tests high_precision_tests.cpp)
add_test(unittests tests)
add_test(hp high_precision_tests)
