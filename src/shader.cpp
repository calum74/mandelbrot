#include "shader.hpp"
#include <cmath>

double fractals::shade(double dx, double dy, double source_x, double source_y, double source_z)
{
    // Copilot created this:
    double distance = std::sqrt(dx * dx + dy * dy);
    double angle = std::atan2(dy, dx);
    double light_angle = std::atan2(source_y - dy, source_x - dx);
    double angle_diff = light_angle - angle;
    double shade = 1.0 - 0.5 * (1.0 + std::cos(angle_diff));
    return shade;
}
