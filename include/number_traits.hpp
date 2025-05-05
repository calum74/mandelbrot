#pragma once

#include <limits>
#include <concepts>
#include <cmath>
#include <utility>

namespace numbers
{
    template<typename T> class number_traits;

    template<std::integral T> struct number_traits<T>
    {
        static constexpr bool is_number = true;
        static constexpr bool is_native = true;
        static constexpr bool is_floating_point = true;
        static constexpr bool is_real = true;
    };

    template<std::floating_point T> struct number_traits<T>
    {
        static constexpr bool is_number = true;
        static constexpr bool is_native = true;
        static constexpr bool is_floating_point = true;
        static constexpr bool is_real = true;
    };

    /*
        Number concepts hierarchy:

        number:
            complex_number: has x.real(), y.real()

            real_number: has to_double(x)
                floating_point:                    
                    native_floating_point
                        use to_double
                    non_native_floating_point: has mantissa_exponent(x);
                fixed_point: has begin(x), end(x), size(x), int_size(x) and value_type
    */

    inline double log(double src);
    inline double to_double(double src) { return src; }
    inline std::pair<double,int> mantissa_exponent(double src) { return {src, 0}; }

    template<typename T>
    concept number = number_traits<T>::is_number;

    template<typename T>
    concept real = number_traits<T>::is_real and requires(T t)
    {
        // May or may not be lossy
        to_double(t);
        log(t);
    };

    template<typename T>
    concept floating_point = number_traits<T>::is_floating_point;

    template<typename T>
    concept native_floating_point = number_traits<T>::is_native;

    template<typename T>
    concept non_native_floating_point = not number_traits<T>::is_native and requires(T t)
    {
        mantissa_exponent(t);
    };

    template<typename T>
    concept fixed_point = real<T> and not floating_point<T> and requires(T t)
    {
        typename number_traits<T>::value_type;
        begin(t);
        end(t);
        size(t);
        int_size(t);
    };
}
