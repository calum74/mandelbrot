#include "RenderingSequence.hpp"

fractals::RenderingSequence::RenderingSequence(int w, int h, int stride)
    : width(w), height(h), initial_stride(stride), stride(0), x(0), y(0) {}

void fractals::RenderingSequence::reset() {
  stride = 0;
  x = 0;
  y = 0;
}

bool fractals::RenderingSequence::next(int &out_x, int &out_y, int &out_stride,
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

bool fractals::RenderingSequence::next0(bool &out_stride_changed) {
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

bool fractals::RenderingSequence::already_done_in_previous_layer() const {
  if (stride == initial_stride)
    return false;

  // Is the current_point already emitted.
  // Look at the least significant bits
  int mask = (stride * 2) - 1;
  return !(x & mask) && !(y & mask);
}
