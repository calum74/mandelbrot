include_directories(../include)
add_library(mandelbrot  experiment.cpp
    mandelbrot.cpp mandeldrop.cpp  view_coords.cpp rendering_sequence.cpp view_parameters.cpp)
