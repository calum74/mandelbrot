#pragma once
#include <iostream>

namespace fractals
{
    // Helper class for representing radii to a higher exponent.

    class radius
    {
    public:
        radius();
        struct from_ln{};

        radius(double ln_r_value, from_ln);
        explicit radius(double value);

        double to_double() const;
        double ln_r() const;
    private:
        double ln_r_value;
    };

    radius operator*(radius a, radius b);
    radius operator/(radius a, radius b);

    bool operator<(radius a, radius b);
    bool operator<=(radius a, radius b);
    bool operator>(radius a, radius b);
    bool operator>=(radius a, radius b);
    bool operator==(radius a, radius b);
    bool operator!=(radius a, radius b);

    // Outputs the radius in engineering format
    std::ostream & operator<<(std::ostream &os, radius r);
}
