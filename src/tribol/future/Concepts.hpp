#ifndef TRIBOL_FUTURE_CONCEPTS_HPP_
#define TRIBOL_FUTURE_CONCEPTS_HPP_

#include "tribol/future/Config.hpp"
#include "tribol/future/Types.hpp"
#include "tribol/future/detail/AutomaticDifferentiation.hpp"

#include <array>
#include <concepts>
#include <limits>
#include <type_traits>

namespace tribol::future {

class ContactDomain;
struct CandidatePairs;

namespace detail {

struct ExecutionProbe {
  TRIBOL_FUTURE_HOST_DEVICE void operator()( int ) const {}
};

}  // namespace detail

template <typename T>
concept ContactConstraint = requires {
  { T::is_nodal } -> std::convertible_to<bool>;
};

template <typename T>
concept ExecutionPolicy = requires( int count, detail::ExecutionProbe probe ) {
  { T::uses_device } -> std::convertible_to<bool>;
  { T::uses_hip } -> std::convertible_to<bool>;
  { T::uses_cuda } -> std::convertible_to<bool>;
  { T::deterministic } -> std::convertible_to<bool>;
  { T::forAll( count, probe ) } -> std::same_as<void>;
  { T::synchronize() } -> std::same_as<void>;
};

template <typename T>
concept SurfaceDiscretization =
    requires( const mfem::ParGridFunction& coordinates, const typename T::Parameters& parameters ) {
      typename T::Parameters;
      { T::is_low_order_refined } -> std::convertible_to<bool>;
      { T::validate( coordinates, parameters ) } -> std::same_as<void>;
    } && ( T::is_low_order_refined || requires {
      { T::max_order } -> std::convertible_to<int>;
    } );

template <typename T>
concept SearchPolicy = std::constructible_from<T, typename T::Parameters> &&
                       requires( T search, const ContactDomain& domain, CandidatePairs& candidates, bool self_contact,
                                 int adjacency_depth ) {
                         typename T::Parameters;
                         {
                           search.findCandidates( domain, candidates, self_contact, adjacency_depth )
                         } -> std::same_as<void>;
                       };

template <typename T>
concept EnforcementLaw =
    std::is_trivially_copyable_v<T> && requires( T law, Real gap, detail::Gradient<Real, 1> differentiable_gap ) {
      { law.potential( gap ) } -> std::convertible_to<Real>;
      { law.pressure( gap ) } -> std::convertible_to<Real>;
      { law.tangent( gap ) } -> std::convertible_to<Real>;
      { law.potential( differentiable_gap ) } -> std::convertible_to<detail::Gradient<Real, 1>>;
      { law.pressure( differentiable_gap ) } -> std::convertible_to<detail::Gradient<Real, 1>>;
      { law.tangent( differentiable_gap ) } -> std::convertible_to<detail::Gradient<Real, 1>>;
    };

template <typename T>
concept CommonPlaneNormalLaw = std::is_trivially_copyable_v<T> && requires( T law, const Real* field ) {
  { T::ties_normal } -> std::convertible_to<bool>;
  { T::ties_tangential } -> std::convertible_to<bool>;
  { T::needs_reference_field } -> std::convertible_to<bool>;
  { T::needs_material_fields } -> std::convertible_to<bool>;
  { law.penaltyStiffness( 0, 1, field, field ) } -> std::convertible_to<Real>;
};

template <typename T>
concept CommonPlaneRateLaw = std::is_trivially_copyable_v<T> && requires( T law, Real stiffness ) {
  { T::needs_velocity } -> std::convertible_to<bool>;
  { T::is_conservative } -> std::convertible_to<bool>;
  { law.rateCoefficient( stiffness ) } -> std::convertible_to<Real>;
};

template <typename T>
concept CommonPlaneTangentialLaw =
    std::is_trivially_copyable_v<T> &&
    requires( T law, std::array<Real, 3> relative_velocity, Real normal_rate, std::array<Real, 3> normal,
              std::array<Real, 3> traction, std::array<detail::Gradient<Real, 1>, 3> differentiable_velocity,
              detail::Gradient<Real, 1> differentiable_rate,
              std::array<detail::Gradient<Real, 1>, 3> differentiable_normal,
              std::array<detail::Gradient<Real, 1>, 3> differentiable_traction ) {
      { T::needs_velocity } -> std::convertible_to<bool>;
      { T::is_conservative } -> std::convertible_to<bool>;
      { law.addTraction( relative_velocity, normal_rate, normal, 3, traction ) } -> std::same_as<void>;
      {
        law.addTraction( differentiable_velocity, differentiable_rate, differentiable_normal, 3,
                         differentiable_traction )
      } -> std::same_as<void>;
    };

template <typename T>
concept MortarBasisPolicy =
    std::is_trivially_copyable_v<T> && requires( T basis, Real coordinate, detail::Gradient<Real, 1> differentiable ) {
      { basis.value( 4, false, false, coordinate, coordinate, coordinate ) } -> std::convertible_to<Real>;
      {
        basis.value( 4, false, false, differentiable, differentiable, differentiable )
      } -> std::convertible_to<detail::Gradient<Real, 1>>;
    };

template <typename T>
concept LinearizationPolicy = requires {
  { T::uses_enzyme } -> std::convertible_to<bool>;
  { T::supports_matrix_free } -> std::convertible_to<bool>;
};

template <typename T>
struct ContactAlgorithmTraits {
  static constexpr MethodCapabilities capabilities = [] {
    if constexpr ( requires {
                     { T::capabilities } -> std::convertible_to<MethodCapabilities>;
                   } ) {
      return MethodCapabilities{ T::capabilities };
    } else {
      return MethodCapabilities{};
    }
  }();

  static constexpr bool accepts_external_pressure = [] {
    if constexpr ( requires {
                     { T::accepts_external_pressure } -> std::convertible_to<bool>;
                   } ) {
      return static_cast<bool>( T::accepts_external_pressure );
    } else {
      return false;
    }
  }();

  static constexpr bool has_nodal_kinematics = [] {
    if constexpr ( requires {
                     { T::has_nodal_kinematics } -> std::convertible_to<bool>;
                   } ) {
      return static_cast<bool>( T::has_nodal_kinematics );
    } else {
      return false;
    }
  }();

  static constexpr int maximum_element_order = [] {
    if constexpr ( requires {
                     { T::maximum_element_order } -> std::convertible_to<int>;
                   } ) {
      return static_cast<int>( T::maximum_element_order );
    } else {
      return std::numeric_limits<int>::max();
    }
  }();
};

template <typename T>
concept ContactAlgorithm = requires( T algorithm, const ContactDomain& domain, const CandidatePairs& candidates,
                                     typename T::InteractionStorage& interactions, const typename T::State& state,
                                     typename T::Kinematics& kinematics, mfem::Vector& residual,
                                     const mfem::Vector& direction, mfem::Vector& derivative ) {
  typename T::Parameters;
  typename T::State;
  typename T::Kinematics;
  typename T::InteractionStorage;
  typename T::Linearization;
  { algorithm.buildInteractions( domain, candidates, interactions ) } -> std::same_as<void>;
  { algorithm.prepare( domain, interactions, kinematics ) } -> std::same_as<void>;
  { algorithm.evaluateKinematics( domain, interactions, state, kinematics ) } -> std::same_as<void>;
  { algorithm.addResidual( domain, interactions, state, kinematics, residual ) } -> std::same_as<Real>;
  {
    algorithm.applyCoordinateDerivative( domain, interactions, state, kinematics, direction, derivative )
  } -> std::same_as<void>;
};

template <typename T>
concept DualContactAlgorithm =
    ContactAlgorithm<T> && ContactAlgorithmTraits<T>::capabilities.needs_dual_space &&
    requires( T algorithm, const ContactDomain& domain, const typename T::InteractionStorage& interactions,
              const typename T::Kinematics& kinematics, const mfem::Vector& direction, mfem::Vector& derivative ) {
      { algorithm.applyDualDerivative( domain, interactions, direction, derivative ) } -> std::same_as<void>;
      { algorithm.applyGapDerivative( domain, interactions, kinematics, direction, derivative ) } -> std::same_as<void>;
    };

template <typename T>
concept NodalKinematicsAlgorithm = ContactAlgorithm<T> && ContactAlgorithmTraits<T>::has_nodal_kinematics &&
                                   requires( const typename T::Kinematics& kinematics ) {
                                     { kinematics.gap } -> std::convertible_to<const mfem::Vector&>;
                                     { kinematics.weighted_gap } -> std::convertible_to<const mfem::Vector&>;
                                     { kinematics.tributary_area } -> std::convertible_to<const mfem::Vector&>;
                                   };

template <typename T>
concept VelocityContactAlgorithm =
    ContactAlgorithm<T> && ContactAlgorithmTraits<T>::capabilities.needs_velocity &&
    requires( T algorithm, const ContactDomain& domain, const typename T::InteractionStorage& interactions,
              const typename T::State& state, const mfem::Vector& direction, mfem::Vector& derivative ) {
      { algorithm.applyVelocityDerivative( domain, interactions, state, direction, derivative ) } -> std::same_as<void>;
    };

template <typename T>
concept OperatorContactAlgorithm =
    ContactAlgorithm<T> && std::constructible_from<T, typename T::Parameters> &&
    ( !ContactAlgorithmTraits<T>::capabilities.needs_dual_space || DualContactAlgorithm<T> ) &&
    ( !ContactAlgorithmTraits<T>::capabilities.needs_velocity || VelocityContactAlgorithm<T> ) &&
    ( !ContactAlgorithmTraits<T>::has_nodal_kinematics || NodalKinematicsAlgorithm<T> );

}  // namespace tribol::future

#endif
