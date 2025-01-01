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

/*
  Compute all points in threads, writing the results to an array of type T.
*/
class async_rendering_sequence {
public:
  async_rendering_sequence(int w, int h, int stride);

  // Blocking call
  void calculate(int threads, std::atomic<bool> &stop);

protected:
  virtual void calculate_point(int x, int y) = 0;
  virtual void layer_complete(int stride) = 0;

  int width, height, stride;
};

template <typename T>
class buffered_rendering_sequence : public async_rendering_sequence {
public:
  buffered_rendering_sequence(int w, int h, int stride)
      : async_rendering_sequence(w, h, stride), output(w * h) {}

protected:
  // !! Unclear whether we need the members to be atomic
  std::vector<std::atomic<T>> output;

  virtual T get_point(int x, int y) = 0;

private:
  void calculate_point(int x, int y) override {
    output[x + y * width] = get_point(x, y);
  }
};

} // namespace fractals
