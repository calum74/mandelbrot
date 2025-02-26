include_directories(../include)
add_library(mandelbrot calculation_pixmap.cpp gradient_stack.cpp
    mandelbrot.cpp mandeldrop.cpp experiment.cpp view_coords.cpp rendering_sequence.cpp shader.cpp view_parameters.cpp)
