cmake_minimum_required(VERSION 3.21)
project(mandelbrot)
set(CMAKE_CXX_STANDARD 17)
enable_testing()

# add_library(mandelbrot
add_executable(tests tests.cpp)
add_test(unittests tests)