*This article is a work in progress.*

# The role of Taylor series in calculating deep zooms of the Mandelbrot set

This article expands on the [Wikipedia article](https://en.wikipedia.org/wiki/Plotting_algorithms_for_the_Mandelbrot_set) on computing the Mandelbrot set. The section we want to focus on is [perturbation theory and series approximation](https://en.wikipedia.org/wiki/Plotting_algorithms_for_the_Mandelbrot_set#Perturbation_theory_and_series_approximation) at the bottom of the article.

Your main enemies when plotting the Mandelbrot set are (a) number of iterations, and (b) arithmetic precision.

When plotting "deep" Mandelbrot sets, we are typically talking about hundreds of thousands, if not millions of iterations per pixel. We also require precision of say $2^{-1000}$ or even higher, and say 1000 bits of precision per number. This means we can't use hardware-based floating point. Both of these problems have the potential to make the naive Mandelbrot calculation infeasibly slow.

## Basic perturbation theory

