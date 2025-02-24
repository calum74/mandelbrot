#pragma once

#include <vector>

namespace fractals {
template <typename T> class pixmap {
public:
  pixmap() : w(0), h(0) {}
  pixmap(int w, int h, T default_value)
      : w(w), h(h), pixels(w * h, default_value) {}

  using size_type = int;

  using iterator = T *;
  iterator begin() { return pixels.data(); }
  iterator end() { return pixels.data() + pixels.size(); }
  T &operator()(int x, int y) { return pixels[x + y * w]; }
  const T &operator()(int x, int y) const { return pixels[x + y * w]; }
  T &operator[](size_type x) { return pixels[x]; }

  size_type size() const { return pixels.size(); }

  size_type width() const { return w; }

  size_type height() const { return h; }

private:
  int w, h;
  std::vector<T> pixels;
};

template <typename T1, typename T2, typename Fn>
void map_pixmap(const pixmap<T1> &src, pixmap<T2> &dest, Fn fn) {
  for (int j = 0; j < src.height(); ++j)
    for (int i = 0; i < src.width(); ++i)
      dest(i, j) = fn(src(i, j));
}

template <typename T, typename Fn>
void map_pixmap(const pixmap<T> &src, pixmap<T> &dest, double dx, double dy,
                double r, Fn fn, T default_value) {
  for (int j = 0; j < dest.height(); ++j)
    for (int i = 0; i < dest.width(); ++i) {
      int i2 = r * i + dx;
      int j2 = r * j + dy;
      auto &to_pixel = dest(i, j);
      if (i2 >= 0 && i2 < dest.width() && j2 >= 0 && j2 < dest.height()) {
        auto &from_pixel = src(i2, j2);
        to_pixel = fn(from_pixel);
      }
      else {
        to_pixel = default_value;
      }
    }
}

template <typename T> struct error_value {
  T value;
  int error;
};

} // namespace fractals
