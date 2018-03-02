#pragma once

#include <boost/range/irange.hpp>

namespace km
{
template <typename T>
auto range(T end)
{
  return boost::irange(static_cast<T>(0), end);
}

template <typename T>
auto range(T start, T end)
{
  return boost::irange(start, end);
}

template <typename T, typename U>
auto range(T start, T end, U step)
{
  return boost::irange(start, end, step);
}
}
