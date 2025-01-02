include_directories(../include)
link_libraries(mandelbrot)

add_library(extra_fractals circle.cpp)

add_executable(tests tests.cpp)
add_executable(high_precision_tests high_precision_tests.cpp)
add_executable(tutorial tutorial.cpp)

add_test(unittests tests)
add_test(hp high_precision_tests)
add_test(tutorial tutorial)
