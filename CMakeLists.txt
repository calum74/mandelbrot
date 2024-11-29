cmake_minimum_required(VERSION 3.21)
project(mandelbrot)
set(CMAKE_CXX_STANDARD 23)
enable_testing()

add_library(mandelbrot center_finder.cpp mandelbrot.cpp registry.cpp)
add_library(extra_fractals circle.cpp)
link_libraries(mandelbrot)

add_executable(tests tests.cpp)
add_test(unittests tests)

