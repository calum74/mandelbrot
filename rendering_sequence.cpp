#include "rendering_sequence.hpp"
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
      if (stride == 1)
        return false;
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

fractals::async_rendering_sequence::async_rendering_sequence(int w, int h,
                                                             int initial_stride)
    : width(w), height(h), stride(initial_stride) {}

fractals::buffered_rendering_sequence::buffered_rendering_sequence(int w, int h,
                                                                   int stride)
    : async_rendering_sequence(w, h, stride), output(w * h) {}

void fractals::async_rendering_sequence::calculate(int threads) {

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
      while (seq.next(x, y, stride, stride_changed)) {
        m.unlock();
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

void fractals::buffered_rendering_sequence::calculate_point(int x, int y) {
  output[x + y * height] = get_point(x, y);
}