// Bits of code that should just be deleted

namespace mandelbrot {

template <typename C> struct orbital {
  C c, z;
  void next() { z = square(z) + c; }
  orbital(C c) : c{c}, z{c} {}
  orbital(C c, C z) : c{c}, z{z} {}
};

template <typename C> orbital<C> mid(orbital<C> o1, orbital<C> o2) {
  return {(o1.c + o2.c) * 0.5, (o1.z + o2.z) * 0.5};
}

// !! Move into mandelbrot.h
template <typename C> int quadrant(C c) { return (c.re > 0) + 2 * (c.im > 0); }

template <typename C> auto scalar_product(C c1, C c2) {
  return c1.re * c2.re + c1.im * c2.im;
}

template <typename C> bool is_right_angle(C c1, C c2) {
  // return false;
  auto sp = scalar_product(c1, c2);
  const double threshold = 0.001;
  return sp * sp < threshold * modulus_squared(c1) * modulus_squared(c2);
}

/*
  A rectangular region.
  The orbits start as a rectangle parallel to the axes, and may become
  scaled and rotated.
  We'll detect when the rectangle gets mapped to something that is no
  longer a rectangle.
*/
struct region {
  using C = mandelbrot::complex<double>;
  // The original region is in `c`
  // The region after `iterations` iterations is `z`
  // C c00, c10, c11, c01, z00, z10, z11, z01;

  // The four orbitals at the corners of the rectangle.
  orbital<C> o00, o01, o11, o10;

  region(orbital<C> a, orbital<C> b, orbital<C> c, orbital<C> d)
      : o00{a}, o10{b}, o11{c}, o01{d} {}

  bool is_rectangle() const {
    // ?? If the region is square, then we could also do right-angles in the
    // center

    return is_right_angle(o10.z - o00.z, o01.z - o00.z) &&
           is_right_angle(o10.z - o00.z, o11.z - o10.z) &&
           is_right_angle(o11.z - o01.z, o11.z - o10.z);
  }

  // We test that the region z is conformal with c
  // This means that z is still a rectangle even after iterating.
  // Is z still a rectangle?
  // We don't test the aspect ratio of the rectangle (maybe we should?)
  bool conformal() const {
    // To what extent is the current region coherent relative to the original
    // region

    // If points are in different quadrants, then we fail
    auto q = quadrant(o00.z);
    if (quadrant(o10.z) != q || quadrant(o11.z) != q || quadrant(o01.z) != q)
      return false;

    return is_rectangle();
  }

  bool fully_escaped() const {
    return escaped(o00.z) && escaped(o10.z) && escaped(o11.z) && escaped(o01.z);
  }

  void iterate() {
    o00.next();
    o01.next();
    o11.next();
    o10.next();
  }
};

} // namespace mandelbrot

class RegionMapperMandelbrot : public NaiveMandelbrot {
  void calculate(fractals::Viewport &view, std::atomic<bool> &stop) override {

    mandelbrot::region r(C{x0, y0}, C{x1, y0}, C{x1, y1}, C{x0, y1});
    regions = 0;
    toplevel = -1;

    calculate_region(view, stop, r, 0, 0, view.width, view.height, 0);
    std::cout << regions << " regions\n";
    std::cout << toplevel << " toplevel iterations\n";
  }

  void set_aspect_ratio(Viewport &vp) override {}

  int max_iterations = 500;
  int toplevel;
  int regions = 0;

  void fill_region(fractals::Viewport &view, int x, int y, int w, int h,
                   int i) {
    ++regions;
    for (int j = 0; j < h; ++j)
      for (int k = 0; k < w; ++k)
        view(x + k, y + j) = fn1(i);
  }

  void calculate_region(fractals::Viewport &view, std::atomic<bool> &stop,
                        mandelbrot::region &r, int x, int y, int w, int h,
                        int i) {

    while (i <= max_iterations && r.conformal() && !r.fully_escaped()) {
      r.iterate();
      ++i;
    }

    if (toplevel <= 0)
      toplevel = i;

    // Compute it fully at each position
    // The region r has given us a starting point
    auto dcx = (r.o10.c - r.o00.c) * (1.0 / w);
    auto dcy = (r.o01.c - r.o00.c) * (1.0 / h);
    auto dzx = (r.o10.z - r.o00.z) * (1.0 / w);
    auto dzy = (r.o01.z - r.o00.z) * (1.0 / h);

    auto max2 = i + 800;

    for (auto rx = 0; rx < w; rx++) {
      for (auto ry = 0; ry < h; ry++) {
        mandelbrot::orbital<C> o{r.o00.c + dcx * (double)rx + dcy * (double)ry,
                                 r.o00.z + dzx * (double)rx + dzy * (double)ry};

        auto i2 = i;
        while (!escaped(o.z) && i2 < max2) {
          o.next();
          ++i2;
        }
        view(x + rx, y + ry) = fn1(i2 == max2 ? 0 : i2);
      }
    }
    return;

    if (r.conformal() && r.fully_escaped()) {
      fill_region(view, x, y, w, h, i);
    } else if (i >= max_iterations) {
      fill_region(view, x, y, w, h, 0);
    } else if (w <= 1 || h <= 1) {
      // Cannot subdivide, so carry on iterating the point
      while (i < max_iterations && !escaped(r.o00.z)) {
        i++;
        r.o00.next();
      }

      fill_region(view, x, y, w, h, i >= max_iterations ? 0 : i);

    } else if (stop) {
      return;
    } else {

      // Split region into four
      auto m0001 = mid(r.o00, r.o01);
      auto m0011 = mid(r.o00, r.o11);
      auto m0010 = mid(r.o00, r.o10);
      auto m1011 = mid(r.o10, r.o11);
      auto m0111 = mid(r.o01, r.o11);

      auto w2 = w / 2;
      auto h2 = h / 2;

      mandelbrot::region r1{r.o00, m0010, m0011, m0001};
      calculate_region(view, stop, r1, x, y, w2, h2, i);

      mandelbrot::region r2{m0010, r.o10, m1011, m0011};
      calculate_region(view, stop, r2, x + w2, y, w - w2, h2, i);

      mandelbrot::region r3{m0011, m1011, r.o11, m0111};
      calculate_region(view, stop, r3, x + w2, y + h2, w - w2, h - h2, i);

      mandelbrot::region r4{m0001, m0011, m0111, r.o01};
      calculate_region(view, stop, r4, x, y + h2, w2, h - h2, i);
    }
  }
};

void report_region2(mandelbrot::region &r, int iterations, int depth) {
  while (r.conformal() && iterations < 100 && !r.fully_escaped()) {
    r.iterate();
    ++iterations;
  }

  std::cout << std::string(depth, ' ') << "Region is conformal to "
            << iterations << " iterations\n";

  if (r.fully_escaped()) {
    std::cout << std::string(depth, ' ') << "Region escaped after "
              << iterations << " iterations\n";
  } else if (depth < 6) {
    auto m0001 = mid(r.o00, r.o01);
    auto m0011 = mid(r.o00, r.o11);
    auto m0010 = mid(r.o00, r.o10);
    auto m1011 = mid(r.o10, r.o11);
    auto m0111 = mid(r.o01, r.o11);

    // Split region into four
    mandelbrot::region r1{r.o00, m0010, m0011, m0001};

    report_region2(r1, iterations, depth + 1);

    mandelbrot::region r2{m0010, r.o10, m1011, m0011};
    report_region2(r2, iterations, depth + 1);

    mandelbrot::region r3{m0011, m1011, r.o11, m0111};
    report_region2(r3, iterations, depth + 1);

    mandelbrot::region r4{m0010, m0011, m0111, r.o01};
    report_region2(r4, iterations, depth + 1);
  }
}

void report_region(C c1, C c2) {
  std::cout << "Region is " << c1.re << " " << c1.im << " " << c2.re << " "
            << c2.im << std::endl;

  mandelbrot::region r(c1, C{c2.re, c1.im}, c2, C{c1.re, c2.im});

  report_region2(r, 0, 0);
}
