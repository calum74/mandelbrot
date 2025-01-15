#pragma once

namespace fractals {
template <typename To, typename From> struct convert_to;

template <typename T> struct convert_to<T, T> {
  static T get(const T &t) { return t; }
};

template <> struct convert_to<double, int> {
  static double get(int i) { return i; }
};

template <typename T1, typename T2> T1 convert(const T2 &x) {
  return convert_to<T1, T2>::get(x);
}

} // namespace fractals