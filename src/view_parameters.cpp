#include "view_parameters.hpp"

#include <iomanip>

std::ostream &fractals::operator<<(std::ostream &os,
                                   const fractals::view_parameters &vp) {

  return os << "Version: 1.00\nX: "
            << std::setprecision(vp.coords.get_precision()) << vp.coords.x
            << "\nY: " << vp.coords.y << "\nR: " << vp.coords.r
            << "\nMax-iterations: " << vp.coords.max_iterations
            << "\nAlgorithm: " << vp.fractal_name
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
      is >> vp.coords.x;
    else if (key == "Y")
      is >> vp.coords.y;
    else if (key == "R")
      is >> vp.coords.r;
    else if (key == "Max-iterations")
      is >> vp.coords.max_iterations;
    else if (key == "Algorithm")
      std::getline(is, vp.fractal_name);
    else if (key == "Colour-scheme")
      is >> vp.colour_seed;
    else if (key == "Colour-gradient")
      is >> vp.colour_gradient;
    else
      break;
  }
  return is;
}
