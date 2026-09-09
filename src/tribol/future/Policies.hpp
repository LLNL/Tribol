#ifndef TRIBOL_FUTURE_POLICIES_HPP_
#define TRIBOL_FUTURE_POLICIES_HPP_

#include "tribol/future/Config.hpp"
#include "tribol/future/detail/AutomaticDifferentiation.hpp"

#include <array>
#include <type_traits>

namespace tribol::future {

struct NodalConstraints {
  static constexpr bool is_nodal = true;
};

struct QuadraturePointConstraints {
  static constexpr bool is_nodal = false;
};

struct QuadraticPenaltyLaw {
  Real stiffness{ 1.0 };

  template <typename T>
  TRIBOL_FUTURE_HOST_DEVICE T potential( const T& gap ) const
  {
    const T active_gap = detail::activeMinimum( gap );
    return T{ 0.5 * stiffness } * active_gap * active_gap;
  }

  template <typename T>
  TRIBOL_FUTURE_HOST_DEVICE T pressure( const T& gap ) const
  {
    return T{ stiffness } * detail::activeMinimum( gap );
  }

  template <typename T>
  TRIBOL_FUTURE_HOST_DEVICE T tangent( const T& gap ) const
  {
    return detail::scalarValue( gap ) < 0.0 ? T{ stiffness } : T{ 0.0 };
  }
};

struct NodalLagrangeMultiplier {
  static constexpr bool is_lagrange_multiplier = true;
};

struct ExternalPressureLaw {
  static constexpr bool is_external = true;
};

struct EnzymeLinearization {
#ifdef TRIBOL_USE_ENZYME
  static constexpr bool available = true;
#else
  static constexpr bool available = false;
#endif
  static constexpr bool uses_enzyme = available;
  static constexpr bool supports_matrix_free = true;
};

struct AnalyticLinearization {
  static constexpr bool uses_enzyme = false;
  static constexpr bool supports_matrix_free = true;
};

struct ConstantKinematicPenalty {
  Real stiffness{ 1.0 };

  static constexpr bool ties_normal = false;
  static constexpr bool ties_tangential = false;
  static constexpr bool needs_reference_field = false;
  static constexpr bool needs_material_fields = false;

  TRIBOL_FUTURE_HOST_DEVICE Real penaltyStiffness( int, int, const Real*, const Real* ) const { return stiffness; }
};

struct MaterialKinematicPenalty {
  Real scale{ 1.0 };

  static constexpr bool ties_normal = false;
  static constexpr bool ties_tangential = false;
  static constexpr bool needs_reference_field = false;
  static constexpr bool needs_material_fields = true;

  TRIBOL_FUTURE_HOST_DEVICE Real penaltyStiffness( int mortar, int nonmortar, const Real* thickness,
                                                   const Real* modulus ) const
  {
    if ( thickness == nullptr || modulus == nullptr ) {
      return scale;
    }
    return scale * 0.5 * ( modulus[mortar] / thickness[mortar] + modulus[nonmortar] / thickness[nonmortar] );
  }
};

struct NoRatePenalty {
  static constexpr bool needs_velocity = false;
  static constexpr bool is_conservative = true;

  TRIBOL_FUTURE_HOST_DEVICE Real rateCoefficient( Real ) const { return 0.0; }
};

struct ConstantRatePenalty {
  Real coefficient{};

  static constexpr bool needs_velocity = true;
  static constexpr bool is_conservative = false;

  TRIBOL_FUTURE_HOST_DEVICE Real rateCoefficient( Real ) const { return coefficient; }
};

struct PercentageRatePenalty {
  Real ratio{};

  static constexpr bool needs_velocity = true;
  static constexpr bool is_conservative = false;

  TRIBOL_FUTURE_HOST_DEVICE Real rateCoefficient( Real normal_stiffness ) const { return ratio * normal_stiffness; }
};

struct FrictionlessTangentialResponse {
  static constexpr bool needs_velocity = false;
  static constexpr bool is_conservative = true;

  template <typename T>
  TRIBOL_FUTURE_HOST_DEVICE void addTraction( const std::array<T, 3>&, const T&, const std::array<T, 3>&, int,
                                              std::array<T, 3>& ) const
  {
  }
};

struct TiedNormalResponse {
  Real stiffness{ 1.0 };

  static constexpr bool ties_normal = true;
  static constexpr bool ties_tangential = false;
  static constexpr bool needs_reference_field = true;
  static constexpr bool needs_material_fields = false;

  TRIBOL_FUTURE_HOST_DEVICE Real penaltyStiffness( int, int, const Real*, const Real* ) const { return stiffness; }
};

struct TiedFullResponse {
  Real stiffness{ 1.0 };

  static constexpr bool ties_normal = true;
  static constexpr bool ties_tangential = true;
  static constexpr bool needs_reference_field = true;
  static constexpr bool needs_material_fields = false;

  TRIBOL_FUTURE_HOST_DEVICE Real penaltyStiffness( int, int, const Real*, const Real* ) const { return stiffness; }
};

struct ViscousTangentialResponse {
  Real damping{};

  static constexpr bool needs_velocity = true;
  static constexpr bool is_conservative = false;

  template <typename T>
  TRIBOL_FUTURE_HOST_DEVICE void addTraction( const std::array<T, 3>& relative_velocity, const T& normal_rate,
                                              const std::array<T, 3>& normal, int dimension,
                                              std::array<T, 3>& traction ) const
  {
    for ( int component = 0; component < dimension; ++component ) {
      traction[component] =
          traction[component] + T{ damping } * ( relative_velocity[component] - normal_rate * normal[component] );
    }
  }
};

struct PrimalMortarBasis {
  template <typename T>
  TRIBOL_FUTURE_HOST_DEVICE T value( int, bool, bool, const T&, const T&, const T& primal_value ) const
  {
    return primal_value;
  }
};

struct DualMortarBasis {
  template <typename T>
  TRIBOL_FUTURE_HOST_DEVICE T value( int nodes, bool high_x, bool high_y, const T& reference_x, const T& reference_y,
                                     const T& primal_value ) const
  {
    if ( nodes == 3 ) {
      return T{ 4.0 } * primal_value - T{ 1.0 };
    }
    const T dual_x = high_x ? T{ 3.0 } * reference_x - T{ 1.0 } : T{ 2.0 } - T{ 3.0 } * reference_x;
    const T dual_y = high_y ? T{ 3.0 } * reference_y - T{ 1.0 } : T{ 2.0 } - T{ 3.0 } * reference_y;
    return dual_x * dual_y;
  }
};

template <typename Enforcement>
inline constexpr bool is_lagrange_multiplier_v = requires { Enforcement::is_lagrange_multiplier; };

template <typename Enforcement>
inline constexpr bool is_external_pressure_v = requires { Enforcement::is_external; };

}  // namespace tribol::future

#endif
