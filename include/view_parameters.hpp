#pragma once
#include <string>
#include <iostream>

namespace fractals {

// View parameters completely determine the image (except for the screen dimensions)
struct view_parameters {
  std::string x, y, r;
  int max_iterations;
  std::string algorithm;
  int colour_seed = 0;
  double colour_gradient = 10;
  std::string title; // User-defined string
};

std::istream &operator>>(std::istream &is, view_parameters &params);
std::ostream &operator<<(std::ostream &is, const view_parameters &params);

bool try_import(std::istream &is, view_parameters &params);
} // namespace fractals
