#include "view_parameters.hpp"

#include <iomanip>
#include <string>

std::ostream &fractals::operator<<(std::ostream &os,
                                   const fractals::view_parameters &vp) {

  return os << "Version: 1.00\nX: "
            << vp.x
            << "\nY: " << vp.y << "\nR: " << vp.r
            << "\nMax-iterations: " << vp.max_iterations
            << "\nAlgorithm: " << vp.algorithm
            << "\nColour-scheme: " << vp.colour_seed
            << "\nColour-gradient: " << vp.colour_gradient << "\n";
}

std::istream &fractals::operator>>(std::istream &is,
                                   fractals::view_parameters &vp) {
  std::string key, value;
  char ch;
  while (is) {
    is >> key;
    is.get();
    key.pop_back();
    if (key == "Version")
      std::getline(is, value);
    else if (key == "X")
      is >> vp.x;
    else if (key == "Y")
      is >> vp.y;
    else if (key == "R")
      is >> vp.r;
    else if (key == "Max-iterations")
      is >> vp.max_iterations;
    else if (key == "Algorithm")
      std::getline(is, vp.algorithm);
    else if (key == "Colour-scheme")
      is >> vp.colour_seed;
    else if (key == "Colour-gradient")
      is >> vp.colour_gradient;
    else
      break;
  }
  return is;
}

bool fractals::try_import(std::istream &is, fractals::view_parameters &vp) {
  return false;
}