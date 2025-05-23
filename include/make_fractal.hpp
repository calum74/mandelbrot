#pragma once
#include "numbers/number_cast.hpp"
#include "fractal_calculation.hpp"
#include "plane.hpp"

namespace fractals {
namespace detail {

template <typename... Ts> class MultiPrecisionFactory;

template <typename T, typename... Ts>
class MultiPrecisionFactory<T, Ts...> : public fractal_calculation_factory {
public:
  MultiPrecisionFactory(std::string name, std::string family)
      : name_{name}, family_{family}, tail{name, family},
        value(std::make_shared<T>()) {}

  view_coords initial_coords() const override { return T::initial_coords(); }
  std::string name() const override { return name_; }

  std::string family() const override { return family_; }

  std::shared_ptr<fractal_calculation>
  create(const view_coords &c, int x, int y,
         std::atomic<bool> &stop) const override {

    if (T::valid_for(c)) {
      value->initialize(c, x, y, stop);
      return value;
    }
    return tail.create(c, x, y, stop);
  }

  bool valid_for(const view_coords &c) const override {
    return T::valid_for(c) || tail.valid_for(c);
  }

private:
  MultiPrecisionFactory<Ts...> tail;

  std::string name_;
  std::string family_;
  std::shared_ptr<fractal_calculation> value;
};

template <typename T>
class MultiPrecisionFactory<T> : public fractal_calculation_factory {
public:
  MultiPrecisionFactory(std::string name, std::string family)
      : name_{name}, family_{family}, value(std::make_shared<T>()) {}

  std::shared_ptr<fractal_calculation>
  create(const view_coords &c, int x, int y,
         std::atomic<bool> &stop) const override {
    value->initialize(c, x, y, stop);
    return value;
  }

  view_coords initial_coords() const override { return T::initial_coords(); }

  std::string name() const override { return name_; }

  std::string family() const override { return family_; }

  bool valid_for(const view_coords &c) const override {
    return T::valid_for(c);
  }

private:
  std::string name_;
  std::string family_;
  std::shared_ptr<fractal_calculation> value;
};

template <typename... Ts> class MultiPrecisionFractal : public fractal {
public:
  MultiPrecisionFractal(std::string name, std::string family)
      : name_(name), family_(family) {}

  std::shared_ptr<fractal_calculation_factory> create() const override {
    return std::make_shared<MultiPrecisionFactory<Ts...>>(name_, family_);
  }

  std::string name() const override { return name_; };
  std::string family() const override { return family_; }

private:
  std::string name_, family_;
};
} // namespace detail

/*
  Creates a fractal which is really a sequence of different fractals which
  can be run at different levels of precision.
*/
template <typename... Ts>
detail::MultiPrecisionFractal<Ts...> make_fractal(const char *name,
                                                  const char *family) {
  return {name, family};
}

template <typename... Ts>
detail::MultiPrecisionFractal<Ts...> make_fractal(const char *name) {
  return {name, name};
}

template <template <int> class Fractal>
detail::MultiPrecisionFractal<
    typename Fractal<50>::type, typename Fractal<128>::type,
    typename Fractal<256>::type, typename Fractal<640>::type,
    typename Fractal<1024>::type, typename Fractal<1596>::type,
    typename Fractal<2048>::type, typename Fractal<4096>::type>
make_fractal(const char *name, const char * family) {
  return {name, family};
}

template <template <int> class Fractal>
auto make_fractal(const char *name) {
  return make_fractal<Fractal>(name, name);
}

} // namespace fractals