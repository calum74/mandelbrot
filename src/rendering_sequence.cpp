#include "rendering_sequence.hpp"
#include <cassert>
#include <future>
#include <mutex>

fractals::rendering_sequence::rendering_sequence(int w, int h, int stride)
    : width(w), height(h), initial_stride(stride), stride(0), x(0), y(0) {}

void fractals::rendering_sequence::reset() {
  stride = 0;
  x = 0;
  y = 0;
}

bool fractals::rendering_sequence::next(int &out_x, int &out_y, int &out_stride,
                                        bool &out_stride_changed) {
  out_stride_changed = false;
  if (stride == -1)
    return false;
  if (stride == 0) {
    stride = initial_stride;
    out_x = x = 0;
    out_y = y = 0;
    out_stride = stride;
    return true;
  }
  while (next0(out_stride_changed)) {
    if (!already_done_in_previous_layer()) {
      out_x = x;
      out_y = y;
      out_stride = stride;
      return true;
    }
  }
  return false;
}

bool fractals::rendering_sequence::next0(bool &out_stride_changed) {
  x += stride;
  if (x >= width) {
    x = 0;
    y += stride;
    if (y >= height) {
      if (stride == 1) {
        stride = -1;
        return false;
      }
      out_stride_changed = true;
      stride = stride / 2;
      x = 0;
      y = 0;
    }
  }

  return true;
}

bool fractals::rendering_sequence::already_done_in_previous_layer() const {
  if (stride == initial_stride)
    return false;

  // Is the current_point already emitted.
  // Look at the least significant bits
  int mask = (stride * 2) - 1;
  return !(x & mask) && !(y & mask);
}

void fractals::rendering_sequence::start_at_stride(int s) {
  stride = s;
  x = 0;
  y = 0;
}

fractals::multi_resolution_sequence::multi_resolution_sequence(int width,
                                                               int height,
                                                               int index) {
  // Calculate the initial block size
  block_size = 1;
  while (block_size < width && block_size < height)
    block_size = block_size * 2;
  int points_in_previous_layer = 0;

  for (; block_size > 0; block_size = block_size / 2) {
    int layer_columns = ((width + block_size - 1) / block_size);
    int layer_rows = ((height + block_size - 1) / block_size);
    int layer_points = layer_rows * layer_columns;

    if (index < layer_points) {
      // We are definitely in this layer, but we don't want to return a
      // point that's already been returned in a previous layer
      bool first_layer = points_in_previous_layer == 0;

      int points_per_even_row = first_layer ? layer_columns : layer_columns / 2;
      int points_per_odd_row = layer_columns;
      int points_per_two_rows = points_per_even_row + points_per_odd_row;

      int layer_index = index - points_in_previous_layer;

      int skip_div = layer_index / points_per_two_rows;
      int skip_mod = layer_index % points_per_two_rows;

      if (skip_mod < points_per_even_row) {
        // We're on an even row
        x = first_layer ? skip_mod * block_size
                        : (1 + 2 * skip_mod) * block_size;
        y = skip_div * block_size * 2;

      } else {
        // We're on an odd row
        x = (skip_mod - points_per_even_row) * block_size;
        y = (skip_div * 2 + 1) * block_size;
      }
      return;
    }

    // Is index in this layer?
    points_in_previous_layer = layer_points;
  }
  // block_size = 0
  x = 0;
  y = 0;
}

fractals::async_rendering_sequence::async_rendering_sequence(int w, int h,
                                                             int initial_stride)
    : width(w), height(h), stride(initial_stride) {}

void fractals::async_rendering_sequence::calculate(int threads,
                                                   std::atomic<bool> &stop) {

  if (threads <= 0)
    threads = std::thread::hardware_concurrency();

  rendering_sequence seq(width, height, stride);
  std::mutex m;

  std::vector<std::future<void>> workers;

  std::vector<int> points_at_stride(stride + 1);

  {
    rendering_sequence tmp_seq(width, height, stride);
    int x, y, s;
    bool c;
    while (tmp_seq.next(x, y, s, c))
      ++points_at_stride[s];
  }

  for (int i = 0; i < threads; ++i) {
    workers.push_back(std::async([&] {
      int x, y, stride;
      bool stride_changed;
      m.lock();
      while (seq.next(x, y, stride, stride_changed) && !stop) {
        m.unlock();
        assert(x >= 0);
        assert(x < width);
        assert(y >= 0);
        assert(y < height);
        calculate_point(x, y);
        m.lock();
        if (--points_at_stride[stride] == 0) {
          layer_complete(stride);
        }
      }
      m.unlock();
    }));
  }

  for (auto &w : workers)
    w.get();
}
