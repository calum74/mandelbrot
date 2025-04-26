include_directories(../include)

add_compile_options(/bigobj)

add_library(mandelbrot 
    calculation_pixmap.cpp
    experiment.cpp 
    gradient_stack.cpp
    mandelbrot.cpp
    mandeldrop.cpp
    radius.cpp
    rendering_sequence.cpp
    shader.cpp
    shader_impl.cpp
    view.cpp
    view_animation.cpp
    view_coords.cpp
    view_parameters.cpp
    view_pixmap.cpp)
