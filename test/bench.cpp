#include "mandelbrot.hpp"
#include <chrono>
#include <sstream>

void benchmark(const char *x, const char *y, const char *r, int iterations,
               int w, int h, const fractals::PointwiseFractal &fractal) {
  auto t1 = std::chrono::high_resolution_clock::now();
  fractals::view_coords coords;
  std::stringstream(x) >> coords.x;
  std::stringstream(y) >> coords.y;
  std::stringstream(r) >> coords.r;
  coords.max_iterations = iterations;

  std::atomic<bool> stop;
  auto calculation = fractal.create(coords, w, h, stop);
  auto t2 = std::chrono::high_resolution_clock::now();
  std::cout << "Setup in " << std::chrono::duration<double>(t2 - t1) << "\n";
  for (int j = 0; j < h; j++)
    for (int i = 0; i < w; i++)
      calculation->calculate(i, j);
  auto t3 = std::chrono::high_resolution_clock::now();
  std::cout << "Total in " << std::chrono::duration<double>(t3 - t1) << "\n";
  std::cout << "Average iterations = " << calculation->average_iterations()
            << "\n";
  std::cout << "Average iterations skipped = " << calculation->average_skipped()
            << "\n";
}

int main(int argc, char **argv) {
  benchmark("-0.5", "0.0", "2.0", 500, 500, 500, mandelbrot_fractal);

  benchmark("-1.40116060900525796", "-0.00000000000113124",
            "0.00000000000000897", 1062122, 100, 100, mandelbrot_fractal);
}
