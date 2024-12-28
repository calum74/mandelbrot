#include "layer.hpp"
#include "fractal.hpp"

fractals::layer::layer(const view_coords &coords, int w, int h,
                       const PointwiseCalculation &calc,
                       std::atomic<bool> &stop, int threads)
    : coords(coords), width{w}, height{h}, points(w * h) {
  for (int j = 0; j < height; ++j) {
    for (int i = 0; i < width; ++i) {
      if (stop)
        return;
      (*this)(i, j) = calc.calculate(i, j);
    }
  }
}
