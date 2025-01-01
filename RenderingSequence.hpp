#pragma once

#include <atomic>
#include <vector>

namespace fractals {
// Visits all points (x,y) in the region of size w*h
class RenderingSequence {
public:
  RenderingSequence() = default; // Invalid state do not use
  RenderingSequence(int w, int h, int stride);

  void reset();

  // Gets the next pixel in the sequence
  // Returns true if succeeded, or false if the sequence is done
  // We can call next() multiple times at the end
  bool next(int &out_x, int &out_y, int &out_stride, bool &out_stride_changed);

private:
  // Visit all points using the stride
  bool next0(bool &out_stride_changed);

  bool already_done_in_previous_layer() const;

  int width, height;
  int x, y; // Current
  int initial_stride;
  int stride;
};

/*
  Compute all points in threads, writing the results to an array of type T.
*/
class async_rendering_sequence {
public:
  // Blocking call
  void calculate_in_threads(int threads);

  virtual double calculate(int x, int y) = 0;
  virtual void layer_completed(int stride) = 0;

protected:
  std::vector<std::atomic<double>> outputs;
};
} // namespace fractals
