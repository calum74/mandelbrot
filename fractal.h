// A simple "C" interface wrapper to a mandelbrot set
#pragma once

union rgb
{
    int i;
    char component[3];
    // red = component[0]
    // green = component[1];
    // blue = component[2];
};

typedef struct Fractal* fractal;

// Receiver for rendered/calculated data.
// pixels could either be an iteration count or encoded RGB.
// pixels contains w*h elements, organised by row. 
typedef void (*render_fn)(void *extra, int x0, int y0, int w, int h, const int * pixels);

// A function from number of iterations to an RGB colour
typedef int (*colourmap_fn)(void * extra, float iterations);

enum FractalOption
{
    fo_max_iteration_count,
    fo_viewport_width,
    fo_viewport_height
};

// Creates a new Mandelbrot set with the default starting position
fractal fract_mandelbrot_new(void);

void fractal_set_option(fractal, enum FractalOption option, int value);
int fractal_get_option(fractal, enum FractalOption option);

void fractal_set_output_size(fractal, int w, int h);
void fractal_set_coordinates(fractal, double x0, double y0, double w, double h);
void fractal_set_colourmap(fractal, colourmap_fn);
// void fractal_set_max_iterations(fractal, int);
// int fractal_get_max_iterations(fractal);

// Zooms in or out about a particular point (cx,cy)
// Zoom in if factor > 1
// Zoom out if factor < 1
void fractal_zoom(fractal, int cx, int cy, double factor);

// Call this to free resources 
void fractal_delete(fractal);

// Calls fn once or multiple times to receive the data.
// `extra` is passed to the render_fn
// Returns when all data rendered.
void fractal_render(fractal, render_fn, void *extra);

// Returns immediately, and calls render_fn on a different thread.
void fractal_render_async(fractal, render_fn, void *extra);

void fractal_render_cancel(fractal);

// Gets the width of the current viewport
double fractal_get_width(fractal);

