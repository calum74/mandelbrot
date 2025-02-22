#pragma once
#include "rendering_sequence.hpp"
#include "pixmap.hpp"

namespace fractals {

    template <typename T>
    struct error_value
    {
        T value;
        int error;
    };

    class calculation_pixmap : public rendering_sequence {
    public:
        std::shared_ptr<fractal_calculation> calculation;
        using value_type = error_value<double>;
        pixmap<value_type> pixels;

        void calculate_point(int x, int y, int w) override;
        std::atomic<double> min_value, max_value;
        std::atomic<std::uint64_t> points_calculated;
    };
}
