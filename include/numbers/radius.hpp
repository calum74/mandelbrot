#pragma once
#include "number_cast.hpp"
#include <iostream>

namespace numbers
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


    template<real From>
    struct number_cast_t<radius, From>
    {
        static radius cast(const From&x) {
            return {log(x), radius::from_ln{}};
      }
    };

    double log(radius r);
    double to_double(radius r);
    std::pair<double, int> mantissa_exponent(radius);

    template<> struct number_traits<radius>
    {
        static constexpr bool is_number = true;
        static constexpr bool is_native = false;
        static constexpr bool is_floating_point = true;
        static constexpr bool is_real = true;
    };
}
