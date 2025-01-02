cmake_minimum_required(VERSION 3.21)
project(mandelbrot)
set(CMAKE_CXX_STANDARD 23)
enable_testing()

add_subdirectory(src)
add_subdirectory(test)

include_directories(include)
link_libraries(mandelbrot)

