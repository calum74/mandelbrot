#pragma once
#include "number_traits.hpp"

namespace numbers
{
  template <typename To, typename From>
  struct number_cast_t;

  template <number To, number From>
  To number_cast(const From & from) { return number_cast_t<To, From>::cast(from); }

  template<native_floating_point To, real From>
  struct number_cast_t<To, From>
  {
    static To cast(const From &from) { return to_double(from); }
  };


} // namespace fractals
