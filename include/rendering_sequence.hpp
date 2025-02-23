#pragma once

#include <atomic>
#include <vector>

namespace fractals {
// Visits all points (x,y) in the region of size w*h
class rendering_sequence {
public:
  rendering_sequence() = default; // Invalid state do not use
  rendering_sequence(int w, int h, int stride);

  void start_at_stride(int s);

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

// Calculates a point in a rendering sequence.
// Organizes each point into a square region
// Renders largest regions first, then progressively renders smaller and smaller
// regions
//
// Input:
// - width: the width of the image to render in pixels
// - height: the height of the image to render in pixels
// - index: the index of the point to render in the range 0..(width*height)-1.
//
// Output:
// - x: the pixel x position of the point to render
// - y: the pixel y position of the point to render
// - block_size: the size of the block to render (may exceeed width and height)
struct multi_resolution_sequence {
  multi_resolution_sequence(int width, int height, int index);
  int x, y, block_size;
};

// Note that asking for a block_size of 0 gives the end of the array
// (width*height)
int get_first_index_in_multi_resolution_sequence(int width, int height,
                                                 int block_size);

/*
  Compute all points in threads, writing the results to an array of type T.
*/
class async_rendering_sequence {
public:
  async_rendering_sequence(int w, int h, int stride);

  // Blocking call
  void calculate(int threads, std::atomic<bool> &stop);

protected:
  virtual void calculate_point(int x, int y, int w) = 0;
  virtual void layer_complete(int stride) = 0;

  int width, height, stride;
};

} // namespace fractals
