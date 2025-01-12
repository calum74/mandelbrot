include_directories(../include)
link_libraries(mandelbrot)

add_library(extra_fractals circle.cpp)

add_executable(tests tests.cpp)
add_executable(arithmetic arithmetic_tests.cpp)
add_executable(tutorial tutorial.cpp)
add_executable(bench bench.cpp)

add_test(unittests tests)
add_test(arithmetic arithmetic)
add_test(tutorial tutorial)
