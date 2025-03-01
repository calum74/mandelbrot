#pragma once
#include "pixmap.hpp"

namespace fractals {

template <typename T> struct error_value {
  T value;
  int error;
};

using view_pixmap = pixmap<error_value<double>>;

// Perform a pixel-by-pixel remapping, interpolate values and don't increase the
// error values
void map_values(const view_pixmap &src, view_pixmap &dest, double dx, double dy,
                double r);

// Perform a pixel-by-pixel remapping, interpolate values and don't increase the
// error values
void interpolate_values(const view_pixmap &src, view_pixmap &dest, double dx,
                        double dy, double r);

void invalidate_values(view_pixmap &pm);

} // namespace fractals