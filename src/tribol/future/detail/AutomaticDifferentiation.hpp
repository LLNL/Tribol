#ifndef TRIBOL_FUTURE_DETAIL_AUTOMATICDIFFERENTIATION_HPP_
#define TRIBOL_FUTURE_DETAIL_AUTOMATICDIFFERENTIATION_HPP_

#include "tribol/future/Config.hpp"

#include <array>
#include <cmath>
#include <type_traits>

namespace tribol::future::detail {

template <typename Value, int NumberOfDerivatives>
struct Gradient {
  Value value{};
  std::array<Value, NumberOfDerivatives> derivative{};

  TRIBOL_FUTURE_HOST_DEVICE Gradient() = default;
  TRIBOL_FUTURE_HOST_DEVICE Gradient( const Value& value_ ) : value( value_ ) {}

  template <typename Scalar>
    requires std::is_arithmetic_v<Scalar>
  TRIBOL_FUTURE_HOST_DEVICE Gradient( Scalar value_ ) : value( Value{ value_ } )
  {
  }

  TRIBOL_FUTURE_HOST_DEVICE static Gradient variable( const Value& value_, int index )
  {
    Gradient result( value_ );
    result.derivative[index] = Value{ 1.0 };
    return result;
  }
};

template <typename T>
TRIBOL_FUTURE_HOST_DEVICE Real scalarValue( const T& value )
{
  if constexpr ( std::is_arithmetic_v<T> ) {
    return static_cast<Real>( value );
  } else {
    return scalarValue( value.value );
  }
}

template <typename T, int N>
TRIBOL_FUTURE_HOST_DEVICE Gradient<T, N> operator+( const Gradient<T, N>& left, const Gradient<T, N>& right )
{
  Gradient<T, N> result( left.value + right.value );
  for ( int i = 0; i < N; ++i ) {
    result.derivative[i] = left.derivative[i] + right.derivative[i];
  }
  return result;
}

template <typename T, int N>
TRIBOL_FUTURE_HOST_DEVICE Gradient<T, N> operator-( const Gradient<T, N>& left, const Gradient<T, N>& right )
{
  Gradient<T, N> result( left.value - right.value );
  for ( int i = 0; i < N; ++i ) {
    result.derivative[i] = left.derivative[i] - right.derivative[i];
  }
  return result;
}

template <typename T, int N>
TRIBOL_FUTURE_HOST_DEVICE Gradient<T, N> operator-( const Gradient<T, N>& value )
{
  Gradient<T, N> result( -value.value );
  for ( int i = 0; i < N; ++i ) {
    result.derivative[i] = -value.derivative[i];
  }
  return result;
}

template <typename T, int N>
TRIBOL_FUTURE_HOST_DEVICE Gradient<T, N> operator*( const Gradient<T, N>& left, const Gradient<T, N>& right )
{
  Gradient<T, N> result( left.value * right.value );
  for ( int i = 0; i < N; ++i ) {
    result.derivative[i] = left.derivative[i] * right.value + left.value * right.derivative[i];
  }
  return result;
}

template <typename T, int N>
TRIBOL_FUTURE_HOST_DEVICE Gradient<T, N> operator/( const Gradient<T, N>& left, const Gradient<T, N>& right )
{
  Gradient<T, N> result( left.value / right.value );
  const T inverse_denominator = T{ 1.0 } / ( right.value * right.value );
  for ( int i = 0; i < N; ++i ) {
    result.derivative[i] =
        ( left.derivative[i] * right.value - left.value * right.derivative[i] ) * inverse_denominator;
  }
  return result;
}

template <typename T, int N>
TRIBOL_FUTURE_HOST_DEVICE Gradient<T, N> sqrt( const Gradient<T, N>& input )
{
  using std::sqrt;
  Gradient<T, N> result( sqrt( input.value ) );
  const T scale = T{ 0.5 } / result.value;
  for ( int i = 0; i < N; ++i ) {
    result.derivative[i] = scale * input.derivative[i];
  }
  return result;
}

template <typename T>
TRIBOL_FUTURE_HOST_DEVICE T activeMinimum( const T& value )
{
  return scalarValue( value ) < 0.0 ? value : T{ 0.0 };
}

template <typename T>
TRIBOL_FUTURE_HOST_DEVICE T clamp( const T& value, Real lower, Real upper )
{
  if ( scalarValue( value ) < lower ) {
    return T{ lower };
  }
  if ( scalarValue( value ) > upper ) {
    return T{ upper };
  }
  return value;
}

}  // namespace tribol::future::detail

#endif
