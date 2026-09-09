#ifndef TRIBOL_FUTURE_RUNTIMEFACTORY_HPP_
#define TRIBOL_FUTURE_RUNTIMEFACTORY_HPP_

#include "tribol/future/MfemContactOperator.hpp"
#include "tribol/future/Search.hpp"
#include "tribol/future/SurfaceDiscretizations.hpp"
#include "tribol/future/algorithms/CommonPlane.hpp"
#include "tribol/future/algorithms/EnergyMortar.hpp"
#include "tribol/future/algorithms/SingleMortar.hpp"

#include <memory>
#include <stdexcept>
#include <variant>

namespace tribol::future {

enum class BuiltInContactMethod
{
  EnergyMortarNodalPenalty,
  EnergyMortarQuadraturePenalty,
  EnergyMortarNodalMultiplier,
  EnergyMortarExternalPressure,
  CommonPlaneFrictionless,
  SingleMortarMultiplier
};

struct BuiltInContactParameters {
  Real penalty_stiffness{ 1.0 };
  int quadrature_points{ 2 };
  Real endpoint_smoothing{};
  Real search_expansion{};
  MfemBridgeParameters bridge{};
};

template <ExecutionPolicy Execution = SequentialExecution, SearchPolicy Search = BvhSearch>
struct BuiltInContactTypes {
  using NodalPenalty = MfemContactOperator<EnergyMortar<NodalConstraints, QuadraticPenaltyLaw>, LowOrderRefinedSurface,
                                           Execution, Search>;
  using QuadraturePenalty = MfemContactOperator<EnergyMortar<QuadraturePointConstraints, QuadraticPenaltyLaw>,
                                                LowOrderRefinedSurface, Execution, Search>;
  using NodalMultiplier = MfemContactOperator<EnergyMortar<NodalConstraints, NodalLagrangeMultiplier>,
                                              LowOrderRefinedSurface, Execution, Search>;
  using ExternalPressure = MfemContactOperator<EnergyMortar<NodalConstraints, ExternalPressureLaw>,
                                               LowOrderRefinedSurface, Execution, Search>;
  using FrictionlessCommonPlane = MfemContactOperator<CommonPlane<>, LowOrderRefinedSurface, Execution, Search>;
  using MultiplierSingleMortar = MfemContactOperator<SingleMortar<>, LowOrderRefinedSurface, Execution, Search>;
  using Variant = std::variant<std::unique_ptr<NodalPenalty>, std::unique_ptr<QuadraturePenalty>,
                               std::unique_ptr<NodalMultiplier>, std::unique_ptr<ExternalPressure>,
                               std::unique_ptr<FrictionlessCommonPlane>, std::unique_ptr<MultiplierSingleMortar>>;
};

template <typename Contact>
typename Contact::Parameters makeContactParameters( const BuiltInContactParameters& parameters )
{
  typename Contact::Parameters result;
  result.search.expansion = parameters.search_expansion;
  result.bridge = parameters.bridge;
  if constexpr ( requires { result.algorithm.quadrature_points; } ) {
    result.algorithm.quadrature_points = parameters.quadrature_points;
  }
  if constexpr ( requires { result.algorithm.endpoint_smoothing; } ) {
    result.algorithm.endpoint_smoothing = parameters.endpoint_smoothing;
  }
  if constexpr ( requires { result.algorithm.enforcement.stiffness; } ) {
    result.algorithm.enforcement.stiffness = parameters.penalty_stiffness;
  }
  if constexpr ( requires { result.algorithm.normal.stiffness; } ) {
    result.algorithm.normal.stiffness = parameters.penalty_stiffness;
  }
  return result;
}

template <ExecutionPolicy Execution = SequentialExecution, SearchPolicy Search = BvhSearch>
typename BuiltInContactTypes<Execution, Search>::Variant makeBuiltInContact(
    BuiltInContactMethod method, const mfem::ParMesh& mesh, const mfem::ParGridFunction& coordinates,
    const PairedSurfaces& surfaces, const BuiltInContactParameters& parameters = {} )
{
  using Types = BuiltInContactTypes<Execution, Search>;
  switch ( method ) {
    case BuiltInContactMethod::EnergyMortarNodalPenalty:
      return std::make_unique<typename Types::NodalPenalty>(
          mesh, coordinates, surfaces, makeContactParameters<typename Types::NodalPenalty>( parameters ) );
    case BuiltInContactMethod::EnergyMortarQuadraturePenalty:
      return std::make_unique<typename Types::QuadraturePenalty>(
          mesh, coordinates, surfaces, makeContactParameters<typename Types::QuadraturePenalty>( parameters ) );
    case BuiltInContactMethod::EnergyMortarNodalMultiplier:
      return std::make_unique<typename Types::NodalMultiplier>(
          mesh, coordinates, surfaces, makeContactParameters<typename Types::NodalMultiplier>( parameters ) );
    case BuiltInContactMethod::EnergyMortarExternalPressure:
      return std::make_unique<typename Types::ExternalPressure>(
          mesh, coordinates, surfaces, makeContactParameters<typename Types::ExternalPressure>( parameters ) );
    case BuiltInContactMethod::CommonPlaneFrictionless:
      return std::make_unique<typename Types::FrictionlessCommonPlane>(
          mesh, coordinates, surfaces, makeContactParameters<typename Types::FrictionlessCommonPlane>( parameters ) );
    case BuiltInContactMethod::SingleMortarMultiplier:
      return std::make_unique<typename Types::MultiplierSingleMortar>(
          mesh, coordinates, surfaces, makeContactParameters<typename Types::MultiplierSingleMortar>( parameters ) );
  }
  throw std::invalid_argument( "Unsupported built-in contact method." );
}

template <ExecutionPolicy Execution = SequentialExecution, SearchPolicy Search = BvhSearch>
typename BuiltInContactTypes<Execution, Search>::Variant makeBuiltInContact(
    BuiltInContactMethod method, const mfem::ParMesh& mesh, const mfem::ParGridFunction& coordinates,
    const SelfContactSurface& surface, const BuiltInContactParameters& parameters = {} )
{
  using Types = BuiltInContactTypes<Execution, Search>;
  if ( method != BuiltInContactMethod::CommonPlaneFrictionless ) {
    throw std::invalid_argument( "Only built-in common-plane contact supports self-contact." );
  }
  return std::make_unique<typename Types::FrictionlessCommonPlane>(
      mesh, coordinates, surface, makeContactParameters<typename Types::FrictionlessCommonPlane>( parameters ) );
}

}  // namespace tribol::future

#endif
