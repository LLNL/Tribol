#ifndef TRIBOL_FUTURE_MFEMCONTACTOPERATOR_HPP_
#define TRIBOL_FUTURE_MFEMCONTACTOPERATOR_HPP_

#include "tribol/future/Concepts.hpp"
#include "tribol/future/Execution.hpp"
#include "tribol/future/MfemBridge.hpp"
#include "tribol/future/Results.hpp"

#include <memory>
#include <stdexcept>
#include <type_traits>

namespace tribol::future {

template <OperatorContactAlgorithm Algorithm, SurfaceDiscretization Surface, ExecutionPolicy Execution,
          SearchPolicy Search>
class MfemContactOperator {
 public:
  using Traits = ContactAlgorithmTraits<Algorithm>;

  struct Parameters {
    typename Algorithm::Parameters algorithm{};
    typename Surface::Parameters surface{};
    typename Search::Parameters search{};
    MfemBridgeParameters bridge{};
  };

  MfemContactOperator( const mfem::ParMesh& mesh, const mfem::ParGridFunction& coordinates,
                       const PairedSurfaces& surfaces, Parameters parameters = {} )
      : bridge_( std::type_identity<Surface>{}, mesh, coordinates, surfaces, parameters.surface,
                 normalizedBridgeParameters( parameters.bridge ) ),
        algorithm_( parameters.algorithm ),
        search_( parameters.search ),
        contact_residual_( bridge_.domain().coordinateSize() ),
        parent_residual_( bridge_.parentCoordinateSpace().GetTrueVSize() ),
        contact_direction_( bridge_.domain().coordinateSize() ),
        contact_derivative_( bridge_.domain().coordinateSize() ),
        mapped_external_pressure_{ external_potential_, external_pressure_, external_tangent_ }
  {
    configureMemory();
  }

  MfemContactOperator( const mfem::ParMesh& mesh, const mfem::ParGridFunction& coordinates,
                       const SelfContactSurface& surface, Parameters parameters = {} )
      : bridge_( std::type_identity<Surface>{}, mesh, coordinates, surface, parameters.surface,
                 normalizedBridgeParameters( parameters.bridge ) ),
        algorithm_( parameters.algorithm ),
        search_( parameters.search ),
        contact_residual_( bridge_.domain().coordinateSize() ),
        parent_residual_( bridge_.parentCoordinateSpace().GetTrueVSize() ),
        contact_direction_( bridge_.domain().coordinateSize() ),
        contact_derivative_( bridge_.domain().coordinateSize() ),
        mapped_external_pressure_{ external_potential_, external_pressure_, external_tangent_ }
  {
    if constexpr ( !Traits::capabilities.supports_self_contact ) {
      throw std::invalid_argument( "The selected contact algorithm does not support self-contact." );
    }
    configureMemory();
  }

  MfemContactOperator( const MfemContactOperator& ) = delete;
  MfemContactOperator& operator=( const MfemContactOperator& ) = delete;
  MfemContactOperator( MfemContactOperator&& ) = delete;
  MfemContactOperator& operator=( MfemContactOperator&& ) = delete;

  void updateInteractions()
  {
    if ( bridge_.maximumSurfaceOrder() > Traits::maximum_element_order ) {
      throw std::invalid_argument(
          "The selected contact algorithm does not support the surface finite-element order; use a compatible "
          "surface discretization." );
    }
    if constexpr ( requires {
                     search_.template findCandidates<Execution>( bridge_.domain(), candidates_, bridge_.selfContact(),
                                                                 bridge_.adjacencyExclusionDepth() );
                   } ) {
      search_.template findCandidates<Execution>( bridge_.domain(), candidates_, bridge_.selfContact(),
                                                  bridge_.adjacencyExclusionDepth() );
    } else {
      search_.findCandidates( bridge_.domain(), candidates_, bridge_.selfContact(), bridge_.adjacencyExclusionDepth() );
    }
    algorithm_.buildInteractions( bridge_.domain(), candidates_, interactions_ );
    algorithm_.prepare( bridge_.domain(), interactions_, kinematics_ );
    bridge_.markInteractionsUpdated();
    prepared_interaction_version_ = bridge_.interactionVersion();
  }

  void updateGeometry( const mfem::ParGridFunction& coordinates ) { bridge_.updateGeometry( coordinates ); }

  void rebuildGeometry( const mfem::ParGridFunction& coordinates )
  {
    bridge_.rebuildGeometry( coordinates );
    resizeWorkspaces();
  }

  NodalKinematicsView evaluateNodalKinematics()
    requires( Traits::has_nodal_kinematics )
  {
    requireInteractions();
    typename Algorithm::State state{};
    evaluateAlgorithmKinematics( state );
    updateTrueNodalKinematics();
    return { true_gap_, true_weighted_gap_, true_tributary_area_ };
  }

  QuadratureDiagnosticsView quadratureDiagnostics() const
  {
    requireInteractions();
    return { kinematics_.pair_offsets, kinematics_.quadrature_points, kinematics_.quadrature_weights,
             kinematics_.quadrature_gaps };
  }

  PenaltyResultView addResidual( mfem::Vector& residual )
    requires( !Traits::capabilities.needs_dual_space && !Traits::accepts_external_pressure )
  {
    typename Algorithm::State state{};
    return addResidualImpl( state, residual );
  }

  PenaltyResultView addResidual( const ExternalPressureData& pressure, mfem::Vector& residual )
    requires( Traits::accepts_external_pressure )
  {
    typename Algorithm::State state{};
    mapExternalPressure( pressure );
    state.external_pressure = &mapped_external_pressure_;
    return addResidualImpl( state, residual );
  }

  PenaltyResultView addResidual( const mfem::Vector& pressure, mfem::Vector& residual )
    requires( Traits::accepts_external_pressure )
  {
    typename Algorithm::State state{};
    mapExternalPressure( pressure );
    state.external_pressure = &mapped_external_pressure_;
    return addResidualImpl( state, residual );
  }

  LagrangeMultiplierResultView addResidual( const mfem::Vector& dual, mfem::Vector& residual )
    requires( Traits::capabilities.needs_dual_space )
  {
    requireInteractions();
    if ( dual.Size() != bridge_.scalarRestriction().Width() ) {
      throw std::invalid_argument( "The multiplier does not match the surface true-dof space." );
    }
    bridge_.scalarRestriction().Mult( dual, contact_dual_ );
    typename Algorithm::State state{};
    state.dual = &contact_dual_;
    evaluateAlgorithmKinematics( state );
    contact_residual_ = 0.0;
    addAlgorithmResidual( state );
    bridge_.coordinateRestriction().MultTranspose( contact_residual_, parent_residual_ );
    addParentResidual( residual );
    updateTrueNodalKinematics();
    return { parent_residual_, true_weighted_gap_ };
  }

  PenaltyResultView addResidual( const typename Algorithm::State& state, mfem::Vector& residual )
    requires( !Traits::capabilities.needs_dual_space )
  {
    return addResidualImpl( state, residual );
  }

  PenaltyLinearization linearize()
    requires( !Traits::capabilities.needs_dual_space && !Traits::accepts_external_pressure )
  {
    typename Algorithm::State state{};
    refreshKinematics( state );
    coordinate_derivative_operator_->reset( state );
    return { coordinate_derivative_operator_.get(), nullptr };
  }

  PenaltyLinearization linearize( const ExternalPressureData& pressure )
    requires( Traits::accepts_external_pressure )
  {
    typename Algorithm::State state{};
    mapExternalPressure( pressure );
    state.external_pressure = &mapped_external_pressure_;
    refreshKinematics( state );
    coordinate_derivative_operator_->reset( state );
    return { coordinate_derivative_operator_.get(), nullptr };
  }

  PenaltyLinearization linearize( const typename Algorithm::State& state )
    requires( !Traits::capabilities.needs_dual_space )
  {
    refreshKinematics( state );
    coordinate_derivative_operator_->reset( state );
    PenaltyLinearization result{ coordinate_derivative_operator_.get(), nullptr };
    if constexpr ( Traits::capabilities.needs_velocity ) {
      velocity_derivative_operator_->reset( state );
      result.dforce_dvelocity = velocity_derivative_operator_.get();
    }
    return result;
  }

  LagrangeMultiplierLinearization linearize( const mfem::Vector& dual )
    requires( Traits::capabilities.needs_dual_space )
  {
    if ( dual.Size() != bridge_.scalarRestriction().Width() ) {
      throw std::invalid_argument( "The multiplier does not match the surface true-dof space." );
    }
    bridge_.scalarRestriction().Mult( dual, contact_dual_ );
    typename Algorithm::State state{};
    state.dual = &contact_dual_;
    refreshKinematics( state );
    coordinate_derivative_operator_->reset( state );
    dual_derivative_operator_->reset();
    gap_derivative_operator_->reset();
    return { coordinate_derivative_operator_.get(), dual_derivative_operator_.get(), gap_derivative_operator_.get() };
  }

  const MfemBridge& bridge() const { return bridge_; }
  const typename Algorithm::InteractionStorage& interactions() const { return interactions_; }
  const typename Algorithm::Kinematics& kinematics() const { return kinematics_; }

 private:
  class CoordinateDerivativeOperator final : public mfem::Operator {
   public:
    explicit CoordinateDerivativeOperator( MfemContactOperator& contact )
        : mfem::Operator( contact.bridge_.parentCoordinateSpace().GetTrueVSize() ), contact_( contact )
    {
    }

    void reset( typename Algorithm::State state )
    {
      state_ = state;
      geometry_version_ = contact_.bridge_.geometryVersion();
      interaction_version_ = contact_.bridge_.interactionVersion();
    }

    void Mult( const mfem::Vector& direction, mfem::Vector& derivative ) const override
    {
      contact_.requireVersion( geometry_version_, interaction_version_ );
      contact_.bridge_.coordinateRestriction().Mult( direction, contact_.contact_direction_ );
      contact_.contact_derivative_ = 0.0;
      contact_.applyAlgorithmCoordinateDerivative( state_ );
      contact_.bridge_.coordinateRestriction().MultTranspose( contact_.contact_derivative_, derivative );
    }

   private:
    MfemContactOperator& contact_;
    typename Algorithm::State state_{};
    GeometryVersion geometry_version_{};
    InteractionVersion interaction_version_{};
  };

  class DualDerivativeOperator final : public mfem::Operator {
   public:
    explicit DualDerivativeOperator( MfemContactOperator& contact )
        : mfem::Operator( contact.bridge_.parentCoordinateSpace().GetTrueVSize(),
                          contact.bridge_.scalarRestriction().Width() ),
          contact_( contact )
    {
    }

    void reset()
    {
      geometry_version_ = contact_.bridge_.geometryVersion();
      interaction_version_ = contact_.bridge_.interactionVersion();
    }

    void Mult( const mfem::Vector& direction, mfem::Vector& derivative ) const override
    {
      if constexpr ( Traits::capabilities.needs_dual_space ) {
        contact_.requireVersion( geometry_version_, interaction_version_ );
        contact_.bridge_.scalarRestriction().Mult( direction, contact_.contact_dual_direction_ );
        contact_.contact_derivative_ = 0.0;
        contact_.applyAlgorithmDualDerivative();
        contact_.bridge_.coordinateRestriction().MultTranspose( contact_.contact_derivative_, derivative );
      }
    }

   private:
    MfemContactOperator& contact_;
    GeometryVersion geometry_version_{};
    InteractionVersion interaction_version_{};
  };

  class GapDerivativeOperator final : public mfem::Operator {
   public:
    explicit GapDerivativeOperator( MfemContactOperator& contact )
        : mfem::Operator( contact.bridge_.scalarRestriction().Width(),
                          contact.bridge_.parentCoordinateSpace().GetTrueVSize() ),
          contact_( contact )
    {
    }

    void reset()
    {
      geometry_version_ = contact_.bridge_.geometryVersion();
      interaction_version_ = contact_.bridge_.interactionVersion();
    }

    void Mult( const mfem::Vector& direction, mfem::Vector& derivative ) const override
    {
      if constexpr ( Traits::capabilities.needs_dual_space ) {
        contact_.requireVersion( geometry_version_, interaction_version_ );
        contact_.bridge_.coordinateRestriction().Mult( direction, contact_.contact_direction_ );
        contact_.applyAlgorithmGapDerivative();
        contact_.bridge_.scalarRestriction().MultTranspose( contact_.contact_gap_derivative_, derivative );
      }
    }

   private:
    MfemContactOperator& contact_;
    GeometryVersion geometry_version_{};
    InteractionVersion interaction_version_{};
  };

  class VelocityDerivativeOperator final : public mfem::Operator {
   public:
    explicit VelocityDerivativeOperator( MfemContactOperator& contact )
        : mfem::Operator( contact.bridge_.parentCoordinateSpace().GetTrueVSize() ), contact_( contact )
    {
    }

    void reset( typename Algorithm::State state )
    {
      state_ = state;
      geometry_version_ = contact_.bridge_.geometryVersion();
      interaction_version_ = contact_.bridge_.interactionVersion();
    }

    void Mult( const mfem::Vector& direction, mfem::Vector& derivative ) const override
    {
      if constexpr ( Traits::capabilities.needs_velocity ) {
        contact_.requireVersion( geometry_version_, interaction_version_ );
        contact_.bridge_.coordinateRestriction().Mult( direction, contact_.contact_direction_ );
        contact_.contact_derivative_ = 0.0;
        contact_.applyAlgorithmVelocityDerivative( state_ );
        contact_.bridge_.coordinateRestriction().MultTranspose( contact_.contact_derivative_, derivative );
      }
    }

   private:
    MfemContactOperator& contact_;
    typename Algorithm::State state_{};
    GeometryVersion geometry_version_{};
    InteractionVersion interaction_version_{};
  };

  static MfemBridgeParameters normalizedBridgeParameters( MfemBridgeParameters parameters )
  {
    parameters.use_device = Execution::uses_device;
    return parameters;
  }

  void configureMemory()
  {
    if constexpr ( Execution::uses_device ) {
      static_assert( Execution::available, "The selected device execution policy is unavailable in this build." );
    }
    resizeWorkspaces();
    coordinate_derivative_operator_ = std::make_unique<CoordinateDerivativeOperator>( *this );
    if constexpr ( Traits::capabilities.needs_dual_space ) {
      dual_derivative_operator_ = std::make_unique<DualDerivativeOperator>( *this );
      gap_derivative_operator_ = std::make_unique<GapDerivativeOperator>( *this );
    }
    if constexpr ( Traits::capabilities.needs_velocity ) {
      velocity_derivative_operator_ = std::make_unique<VelocityDerivativeOperator>( *this );
    }
  }

  void resizeWorkspaces()
  {
    contact_residual_.SetSize( bridge_.domain().coordinateSize() );
    contact_direction_.SetSize( bridge_.domain().coordinateSize() );
    contact_derivative_.SetSize( bridge_.domain().coordinateSize() );
    parent_residual_.SetSize( bridge_.parentCoordinateSpace().GetTrueVSize() );
    if constexpr ( Traits::capabilities.needs_dual_space ) {
      contact_dual_.SetSize( bridge_.scalarRestriction().Height() );
      contact_dual_direction_.SetSize( bridge_.scalarRestriction().Height() );
      contact_gap_derivative_.SetSize( bridge_.scalarRestriction().Height() );
      true_gap_derivative_.SetSize( bridge_.scalarRestriction().Width() );
    }
    if constexpr ( Traits::has_nodal_kinematics ) {
      true_gap_.SetSize( bridge_.scalarRestriction().Width() );
      true_weighted_gap_.SetSize( bridge_.scalarRestriction().Width() );
      true_tributary_area_.SetSize( bridge_.scalarRestriction().Width() );
    }
    if constexpr ( Traits::accepts_external_pressure ) {
      external_potential_.SetSize( bridge_.scalarRestriction().Height() );
      external_pressure_.SetSize( bridge_.scalarRestriction().Height() );
      external_tangent_.SetSize( bridge_.scalarRestriction().Height() );
    }
    contact_residual_.UseDevice( Execution::uses_device );
    contact_direction_.UseDevice( Execution::uses_device );
    contact_derivative_.UseDevice( Execution::uses_device );
    parent_residual_.UseDevice( Execution::uses_device );
    if constexpr ( Traits::capabilities.needs_dual_space ) {
      contact_dual_.UseDevice( Execution::uses_device );
      contact_dual_direction_.UseDevice( Execution::uses_device );
      contact_gap_derivative_.UseDevice( Execution::uses_device );
      true_gap_derivative_.UseDevice( Execution::uses_device );
    }
    if constexpr ( Traits::has_nodal_kinematics ) {
      true_gap_.UseDevice( Execution::uses_device );
      true_weighted_gap_.UseDevice( Execution::uses_device );
      true_tributary_area_.UseDevice( Execution::uses_device );
    }
    if constexpr ( Traits::accepts_external_pressure ) {
      external_potential_.UseDevice( Execution::uses_device );
      external_pressure_.UseDevice( Execution::uses_device );
      external_tangent_.UseDevice( Execution::uses_device );
    }
  }

  void requireInteractions() const
  {
    if ( !bridge_.topologyIsCurrent() ) {
      throw std::logic_error(
          "The parent mesh or finite-element space changed; construct a new contact operator after AMR or "
          "repartitioning." );
    }
    if ( !bridge_.interactionsValid() || prepared_interaction_version_ != bridge_.interactionVersion() ) {
      throw std::logic_error(
          "Contact interactions are unavailable; call updateInteractions() after construction or rebuildGeometry()." );
    }
  }

  void requireVersion( GeometryVersion geometry_version, InteractionVersion interaction_version ) const
  {
    requireInteractions();
    if ( geometry_version != bridge_.geometryVersion() || interaction_version != bridge_.interactionVersion() ) {
      throw std::logic_error( "A contact linearization was invalidated by a geometry or interaction update." );
    }
  }

  void refreshKinematics( const typename Algorithm::State& state )
  {
    requireInteractions();
    evaluateAlgorithmKinematics( state );
  }

  void evaluateAlgorithmKinematics( const typename Algorithm::State& state )
  {
    if constexpr ( requires {
                     algorithm_.template evaluateKinematics<Execution>( bridge_.domain(), interactions_, state,
                                                                        kinematics_ );
                   } ) {
      algorithm_.template evaluateKinematics<Execution>( bridge_.domain(), interactions_, state, kinematics_ );
    } else {
      algorithm_.evaluateKinematics( bridge_.domain(), interactions_, state, kinematics_ );
    }
  }

  Real addAlgorithmResidual( const typename Algorithm::State& state )
  {
    if constexpr ( requires {
                     algorithm_.template addResidual<Execution>( bridge_.domain(), interactions_, state, kinematics_,
                                                                 contact_residual_ );
                   } ) {
      return algorithm_.template addResidual<Execution>( bridge_.domain(), interactions_, state, kinematics_,
                                                         contact_residual_ );
    } else {
      return algorithm_.addResidual( bridge_.domain(), interactions_, state, kinematics_, contact_residual_ );
    }
  }

  void applyAlgorithmCoordinateDerivative( const typename Algorithm::State& state )
  {
    if constexpr ( requires {
                     algorithm_.template applyCoordinateDerivative<Execution>(
                         bridge_.domain(), interactions_, state, kinematics_, contact_direction_, contact_derivative_ );
                   } ) {
      algorithm_.template applyCoordinateDerivative<Execution>( bridge_.domain(), interactions_, state, kinematics_,
                                                                contact_direction_, contact_derivative_ );
    } else {
      algorithm_.applyCoordinateDerivative( bridge_.domain(), interactions_, state, kinematics_, contact_direction_,
                                            contact_derivative_ );
    }
  }

  void applyAlgorithmDualDerivative()
    requires( Traits::capabilities.needs_dual_space )
  {
    if constexpr ( requires {
                     algorithm_.template applyDualDerivative<Execution>( bridge_.domain(), interactions_,
                                                                         contact_dual_direction_, contact_derivative_ );
                   } ) {
      algorithm_.template applyDualDerivative<Execution>( bridge_.domain(), interactions_, contact_dual_direction_,
                                                          contact_derivative_ );
    } else {
      algorithm_.applyDualDerivative( bridge_.domain(), interactions_, contact_dual_direction_, contact_derivative_ );
    }
  }

  void applyAlgorithmGapDerivative()
    requires( Traits::capabilities.needs_dual_space )
  {
    if constexpr ( requires {
                     algorithm_.template applyGapDerivative<Execution>( bridge_.domain(), interactions_, kinematics_,
                                                                        contact_direction_, contact_gap_derivative_ );
                   } ) {
      algorithm_.template applyGapDerivative<Execution>( bridge_.domain(), interactions_, kinematics_,
                                                         contact_direction_, contact_gap_derivative_ );
    } else {
      algorithm_.applyGapDerivative( bridge_.domain(), interactions_, kinematics_, contact_direction_,
                                     contact_gap_derivative_ );
    }
  }

  void applyAlgorithmVelocityDerivative( const typename Algorithm::State& state )
    requires( Traits::capabilities.needs_velocity )
  {
    if constexpr ( requires {
                     algorithm_.template applyVelocityDerivative<Execution>( bridge_.domain(), interactions_, state,
                                                                             contact_direction_, contact_derivative_ );
                   } ) {
      algorithm_.template applyVelocityDerivative<Execution>( bridge_.domain(), interactions_, state,
                                                              contact_direction_, contact_derivative_ );
    } else {
      algorithm_.applyVelocityDerivative( bridge_.domain(), interactions_, state, contact_direction_,
                                          contact_derivative_ );
    }
  }

  void updateTrueNodalKinematics()
    requires( Traits::has_nodal_kinematics )
  {
    bridge_.scalarRestriction().MultTranspose( kinematics_.weighted_gap, true_weighted_gap_ );
    bridge_.scalarRestriction().MultTranspose( kinematics_.tributary_area, true_tributary_area_ );
    const bool use_device = Execution::uses_device;
    const Real* weighted_gap = true_weighted_gap_.Read( use_device );
    const Real* area = true_tributary_area_.Read( use_device );
    Real* gap = true_gap_.Write( use_device );
    Execution::forAll( true_gap_.Size(), [=] MFEM_HOST_DEVICE( int dof ) {
      gap[dof] = area[dof] > 0.0 ? weighted_gap[dof] / area[dof] : 0.0;
    } );
  }

  void mapExternalPressure( const ExternalPressureData& pressure )
    requires( Traits::accepts_external_pressure )
  {
    const int expected_size = bridge_.scalarRestriction().Width();
    if ( pressure.potential_density.Size() != expected_size || pressure.pressure.Size() != expected_size ||
         pressure.pressure_tangent.Size() != expected_size ) {
      throw std::invalid_argument( "External pressure data must match the surface true-dof space." );
    }
    bridge_.scalarRestriction().Mult( pressure.potential_density, external_potential_ );
    bridge_.scalarRestriction().Mult( pressure.pressure, external_pressure_ );
    bridge_.scalarRestriction().Mult( pressure.pressure_tangent, external_tangent_ );
  }

  void mapExternalPressure( const mfem::Vector& pressure )
    requires( Traits::accepts_external_pressure )
  {
    if ( pressure.Size() != bridge_.scalarRestriction().Width() ) {
      throw std::invalid_argument( "External pressure must match the surface true-dof space." );
    }
    external_potential_ = 0.0;
    external_tangent_ = 0.0;
    bridge_.scalarRestriction().Mult( pressure, external_pressure_ );
  }

  PenaltyResultView addResidualImpl( const typename Algorithm::State& state, mfem::Vector& residual )
  {
    refreshKinematics( state );
    contact_residual_ = 0.0;
    const Real energy = addAlgorithmResidual( state );
    bridge_.coordinateRestriction().MultTranspose( contact_residual_, parent_residual_ );
    addParentResidual( residual );
    return { parent_residual_, energy };
  }

  void addParentResidual( mfem::Vector& residual ) const
  {
    if ( residual.Size() != parent_residual_.Size() ) {
      throw std::invalid_argument( "The host residual does not match the parent coordinate true-dof space." );
    }
    residual += parent_residual_;
  }

  MfemBridge bridge_;
  Algorithm algorithm_;
  Search search_;
  CandidatePairs candidates_;
  typename Algorithm::InteractionStorage interactions_;
  typename Algorithm::Kinematics kinematics_;
  mfem::Vector contact_residual_;
  mfem::Vector parent_residual_;
  mfem::Vector contact_direction_;
  mfem::Vector contact_derivative_;
  mfem::Vector contact_dual_;
  mfem::Vector contact_dual_direction_;
  mfem::Vector contact_gap_derivative_;
  mfem::Vector true_gap_;
  mfem::Vector true_weighted_gap_;
  mfem::Vector true_tributary_area_;
  mfem::Vector true_gap_derivative_;
  mfem::Vector external_potential_;
  mfem::Vector external_pressure_;
  mfem::Vector external_tangent_;
  ExternalPressureData mapped_external_pressure_;
  InteractionVersion prepared_interaction_version_;
  std::unique_ptr<CoordinateDerivativeOperator> coordinate_derivative_operator_;
  std::unique_ptr<DualDerivativeOperator> dual_derivative_operator_;
  std::unique_ptr<GapDerivativeOperator> gap_derivative_operator_;
  std::unique_ptr<VelocityDerivativeOperator> velocity_derivative_operator_;
};

}  // namespace tribol::future

#endif
