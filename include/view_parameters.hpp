#pragma once
#include <string>
#include <iostream>
#include "shader_parameters.hpp"

namespace fractals {

// View parameters completely determine the image (except for the screen dimensions)
struct view_parameters {
  std::string x, y, r;
  int max_iterations;
  std::string algorithm;
  std::string title; // User-defined string

  shader_parameters shader;
};

std::istream &operator>>(std::istream &is, view_parameters &params);
std::ostream &operator<<(std::ostream &is, const view_parameters &params);

bool try_import(std::istream &is, view_parameters &params);
} // namespace fractals
