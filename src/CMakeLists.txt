include_directories(../include)
add_library(mandelbrot 
    calculation_pixmap.cpp
    experiment.cpp 
    gradient_stack.cpp
    mandelbrot.cpp
    mandeldrop.cpp
    rendering_sequence.cpp
    shader.cpp
    view.cpp
    view_animation.cpp
    view_coords.cpp
    view_parameters.cpp
    view_pixmap.cpp)
