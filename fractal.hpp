#include <cstddef>

/*
Basic infrastructure for rendering fractals.
*/

namespace fractals
{
    struct pixel
    {
        std::uint32_t colour;
    };

    class point
    {
        std::int64_t mantissa, exponent;
    };

    class renderer
    {
    };

    struct region
    {
        point origin, size;
    };

    class viewport
    {
        int width, height;
    };

    class viewer
    {
        virtual void on_update(int x0, int x0, int w, int h, const pixel * output);
    };

    class fractal
    {
    public:
        virtual ~fractal() = default;
        virtual void calculate(viewer);
    };
};
