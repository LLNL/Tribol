#ifndef TRIBOL_FUTURE_DETAIL_ELEMENTGEOMETRY_HPP_
#define TRIBOL_FUTURE_DETAIL_ELEMENTGEOMETRY_HPP_

#include "tribol/future/Config.hpp"
#include "tribol/future/detail/AutomaticDifferentiation.hpp"

#include <array>

namespace tribol::future::detail {

template <typename T, std::size_t Dimension>
using Vector = std::array<T, Dimension>;

template <typename T, std::size_t Dimension>
TRIBOL_FUTURE_HOST_DEVICE Vector<T, Dimension> operator+( const Vector<T, Dimension>& left,
                                                          const Vector<T, Dimension>& right )
{
  Vector<T, Dimension> result{};
  for ( std::size_t component = 0; component < Dimension; ++component ) {
    result[component] = left[component] + right[component];
  }
  return result;
}

template <typename T, std::size_t Dimension>
TRIBOL_FUTURE_HOST_DEVICE Vector<T, Dimension> operator-( const Vector<T, Dimension>& left,
                                                          const Vector<T, Dimension>& right )
{
  Vector<T, Dimension> result{};
  for ( std::size_t component = 0; component < Dimension; ++component ) {
    result[component] = left[component] - right[component];
  }
  return result;
}

template <typename T, std::size_t Dimension>
TRIBOL_FUTURE_HOST_DEVICE Vector<T, Dimension> operator*( const T& scalar, const Vector<T, Dimension>& vector )
{
  Vector<T, Dimension> result{};
  for ( std::size_t component = 0; component < Dimension; ++component ) {
    result[component] = scalar * vector[component];
  }
  return result;
}

template <typename T, std::size_t Dimension>
TRIBOL_FUTURE_HOST_DEVICE T dot( const Vector<T, Dimension>& left, const Vector<T, Dimension>& right )
{
  T result{};
  for ( std::size_t component = 0; component < Dimension; ++component ) {
    result = result + left[component] * right[component];
  }
  return result;
}

template <typename T, std::size_t Dimension>
TRIBOL_FUTURE_HOST_DEVICE T norm( const Vector<T, Dimension>& vector )
{
  using std::sqrt;
  return sqrt( dot( vector, vector ) );
}

template <typename T>
TRIBOL_FUTURE_HOST_DEVICE T cross2d( const Vector<T, 2>& left, const Vector<T, 2>& right )
{
  return left[0] * right[1] - left[1] * right[0];
}

template <typename T>
TRIBOL_FUTURE_HOST_DEVICE Vector<T, 2> edgeNormal( const Vector<T, 2>& first, const Vector<T, 2>& second )
{
  const auto tangent = second - first;
  const auto length = norm( tangent );
  return { tangent[1] / length, -tangent[0] / length };
}

template <typename T>
TRIBOL_FUTURE_HOST_DEVICE T projectAlongNormal( const Vector<T, 2>& point, const Vector<T, 2>& edge_origin,
                                                const Vector<T, 2>& edge_tangent,
                                                const Vector<T, 2>& projection_normal )
{
  return cross2d( point - edge_origin, projection_normal ) / cross2d( edge_tangent, projection_normal );
}

}  // namespace tribol::future::detail

#endif
