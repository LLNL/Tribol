#ifndef TRIBOL_FUTURE_ALGORITHMS_ENERGYMORTAR_HPP_
#define TRIBOL_FUTURE_ALGORITHMS_ENERGYMORTAR_HPP_

#include "tribol/future/Concepts.hpp"
#include "tribol/future/ContactDomain.hpp"
#include "tribol/future/Execution.hpp"
#include "tribol/future/InteractionBatches.hpp"
#include "tribol/future/Policies.hpp"
#include "tribol/future/Results.hpp"
#include "tribol/future/detail/AutomaticDifferentiation.hpp"
#include "tribol/future/detail/ElementGeometry.hpp"
#include "tribol/future/detail/Enzyme.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace tribol::future {

template <ContactConstraint ConstraintPolicy, typename EnforcementPolicy,
          LinearizationPolicy LinearizationPolicy = EnzymeLinearization>
  requires( EnforcementLaw<EnforcementPolicy> || is_lagrange_multiplier_v<EnforcementPolicy> ||
            is_external_pressure_v<EnforcementPolicy> )
class EnergyMortar {
 public:
  static constexpr int max_nodes_per_segment = native_high_order_max_order + 1;
  static constexpr int max_local_coordinates = 4 * max_nodes_per_segment;
  static constexpr int maximum_element_order = native_high_order_max_order;

  struct Parameters {
    EnforcementPolicy enforcement{};
    int quadrature_points{ 2 };
    int projection_iterations{ 6 };
    Real endpoint_smoothing{};
    Real minimum_overlap{ 1.0e-14 };
  };

  struct State {
    const mfem::Vector* dual{};
    const ExternalPressureData* external_pressure{};
  };

  struct Kinematics {
    mfem::Vector gap;
    mfem::Vector weighted_gap;
    mfem::Vector tributary_area;
    mfem::Vector pressure;
    mfem::Vector pressure_tangent;
    mfem::Vector potential_density;
    mfem::Array<int> pair_offsets;
    mfem::Vector quadrature_points;
    mfem::Vector quadrature_weights;
    mfem::Vector quadrature_gaps;
    mutable mfem::Vector directional_weighted_gap;
    mutable mfem::Vector directional_area;
    mutable mfem::Vector scalar_workspace;
  };

  struct InteractionStorage {
    mfem::Array<ElementPair> pairs;
    mfem::Array<InteractionBatch> batches;
    InteractionColoring coloring;
  };

  struct Linearization {
    using policy_type = LinearizationPolicy;
  };

  static constexpr bool accepts_external_pressure = is_external_pressure_v<EnforcementPolicy>;
  static constexpr bool has_nodal_kinematics = ConstraintPolicy::is_nodal;
  static constexpr MethodCapabilities capabilities{
      .needs_dual_space = is_lagrange_multiplier_v<EnforcementPolicy>,
      .needs_velocity = false,
      .needs_reference_field = false,
      .needs_material_fields = false,
      .has_history = false,
      .has_energy = !is_lagrange_multiplier_v<EnforcementPolicy>,
      .supports_self_contact = false,
      .supports_assembled_jacobian = true,
  };

  explicit EnergyMortar( Parameters parameters = {} ) : parameters_( parameters )
  {
    if ( parameters_.quadrature_points < 1 || parameters_.quadrature_points > 3 ) {
      throw std::invalid_argument( "EnergyMortar supports one to three Gauss points." );
    }
    if ( parameters_.projection_iterations < 1 ) {
      throw std::invalid_argument( "EnergyMortar requires at least one fixed projection iteration." );
    }
    if ( parameters_.endpoint_smoothing < 0.0 || parameters_.endpoint_smoothing >= 0.5 ) {
      throw std::invalid_argument( "EnergyMortar endpoint smoothing must be in [0, 0.5)." );
    }
    if constexpr ( !is_lagrange_multiplier_v<EnforcementPolicy> && !is_external_pressure_v<EnforcementPolicy> ) {
      if ( parameters_.enforcement.stiffness <= 0.0 ) {
        throw std::invalid_argument( "EnergyMortar penalty stiffness must be positive." );
      }
    }
    static_assert( ConstraintPolicy::is_nodal || !is_lagrange_multiplier_v<EnforcementPolicy>,
                   "Quadrature-point Lagrange multipliers are not supported." );
    static_assert( ConstraintPolicy::is_nodal || !is_external_pressure_v<EnforcementPolicy>,
                   "External pressure is supported only by nodal EnergyMortar." );
  }

  void buildInteractions( const ContactDomain& domain, const CandidatePairs& candidates,
                          InteractionStorage& interactions ) const
  {
    if ( domain.dimension() != 2 ) {
      throw std::invalid_argument( "EnergyMortar currently supports only two-dimensional surface contact." );
    }
    const auto view = domain.view( false );
    std::vector<ElementPair> active;
    active.reserve( candidates.pairs.Size() );
    for ( int i = 0; i < candidates.pairs.Size(); ++i ) {
      const auto pair = candidates.pairs[i];
      if ( view.topology( pair.mortar_element ) != ElementTopology::Segment ||
           view.topology( pair.nonmortar_element ) != ElementTopology::Segment ) {
        continue;
      }
      if ( view.numberOfElementDofs( pair.mortar_element ) > max_nodes_per_segment ||
           view.numberOfElementDofs( pair.nonmortar_element ) > max_nodes_per_segment ) {
        throw std::invalid_argument( "EnergyMortar segment order exceeds TRIBOL_FUTURE_NATIVE_HO_MAX_ORDER." );
      }
      if ( overlapLength( view, pair ) > parameters_.minimum_overlap ) {
        active.push_back( pair );
      }
    }
    buildInteractionBatches( domain, std::move( active ), interactions.pairs, interactions.batches );
    buildInteractionColoring( domain, interactions.pairs, interactions.coloring );
  }

  void prepare( const ContactDomain& domain, const InteractionStorage& interactions, Kinematics& kinematics ) const
  {
    const int number_of_dofs = domain.numberOfDofs();
    kinematics.gap.SetSize( number_of_dofs );
    kinematics.weighted_gap.SetSize( number_of_dofs );
    kinematics.tributary_area.SetSize( number_of_dofs );
    kinematics.pressure.SetSize( number_of_dofs );
    kinematics.pressure_tangent.SetSize( number_of_dofs );
    kinematics.potential_density.SetSize( number_of_dofs );
    kinematics.directional_weighted_gap.SetSize( number_of_dofs );
    kinematics.directional_area.SetSize( number_of_dofs );
    kinematics.scalar_workspace.SetSize( interactions.pairs.Size() + 1 );
    kinematics.pair_offsets.SetSize( interactions.pairs.Size() + 1 );
    for ( int pair = 0; pair <= interactions.pairs.Size(); ++pair ) {
      kinematics.pair_offsets[pair] = pair * parameters_.quadrature_points;
    }
    const int number_of_points = interactions.pairs.Size() * parameters_.quadrature_points;
    kinematics.quadrature_points.SetSize( number_of_points * 4 );
    kinematics.quadrature_weights.SetSize( number_of_points );
    kinematics.quadrature_gaps.SetSize( number_of_points );
  }

  template <ExecutionPolicy Execution>
  void evaluateKinematics( const ContactDomain& domain, const InteractionStorage& interactions, const State& state,
                           Kinematics& kinematics ) const
  {
    const bool use_device = Execution::uses_device;
    setMemoryMode( kinematics, interactions, use_device );
    kinematics.gap = 0.0;
    kinematics.weighted_gap = 0.0;
    kinematics.tributary_area = 0.0;
    kinematics.pressure = 0.0;
    kinematics.pressure_tangent = 0.0;
    kinematics.potential_density = 0.0;
    kinematics.quadrature_points = 0.0;
    kinematics.quadrature_weights = 0.0;
    kinematics.quadrature_gaps = 0.0;

    const auto view = domain.view( use_device );
    const auto* pairs = interactions.pairs.Read( use_device );
    auto* weighted_gap = kinematics.weighted_gap.ReadWrite( use_device );
    auto* area = kinematics.tributary_area.ReadWrite( use_device );
    auto* quadrature_points = kinematics.quadrature_points.Write( use_device );
    auto* quadrature_weights = kinematics.quadrature_weights.Write( use_device );
    auto* quadrature_gaps = kinematics.quadrature_gaps.Write( use_device );
    const auto parameters = parameters_;

    forEachInteraction<Execution>( interactions.batches, interactions.coloring, [=] MFEM_HOST_DEVICE( int pair_index ) {
      const auto values = pairValues<Real>( view, pairs[pair_index], parameters );
      const int mortar_element = pairs[pair_index].mortar_element;
      const int mortar_nodes = view.numberOfElementDofs( mortar_element );
      for ( int node = 0; node < mortar_nodes; ++node ) {
        const int dof = view.elementDof( mortar_element, node );
        AtomicAdd( weighted_gap[dof], values.weighted_gap[node] );
        AtomicAdd( area[dof], values.area[node] );
      }
      for ( int q = 0; q < parameters.quadrature_points; ++q ) {
        const int point = pair_index * parameters.quadrature_points + q;
        quadrature_gaps[point] = values.gaps[q];
        quadrature_weights[point] = values.weights[q];
        for ( int component = 0; component < 2; ++component ) {
          quadrature_points[( 2 * point ) * 2 + component] = values.mortar_points[q][component];
          quadrature_points[( 2 * point + 1 ) * 2 + component] = values.nonmortar_points[q][component];
        }
      }
    } );

    auto* gap = kinematics.gap.ReadWrite( use_device );
    const int number_of_dofs = domain.numberOfDofs();
    Execution::forAll( number_of_dofs, [=] MFEM_HOST_DEVICE( int dof ) {
      gap[dof] = area[dof] > 0.0 ? weighted_gap[dof] / area[dof] : 0.0;
    } );

    if constexpr ( is_lagrange_multiplier_v<EnforcementPolicy> ) {
      if ( state.dual != nullptr ) {
        kinematics.pressure = *state.dual;
      }
    } else if constexpr ( is_external_pressure_v<EnforcementPolicy> ) {
      if ( state.external_pressure != nullptr ) {
        kinematics.pressure = state.external_pressure->pressure;
        kinematics.pressure_tangent = state.external_pressure->pressure_tangent;
        kinematics.potential_density = state.external_pressure->potential_density;
      }
    } else {
      const auto law = parameters_.enforcement;
      auto* pressure = kinematics.pressure.ReadWrite( use_device );
      auto* tangent = kinematics.pressure_tangent.ReadWrite( use_device );
      auto* potential = kinematics.potential_density.ReadWrite( use_device );
      Execution::forAll( number_of_dofs, [=] MFEM_HOST_DEVICE( int dof ) {
        pressure[dof] = law.pressure( gap[dof] );
        tangent[dof] = law.tangent( gap[dof] );
        potential[dof] = law.potential( gap[dof] );
      } );
    }
  }

  void evaluateKinematics( const ContactDomain& domain, const InteractionStorage& interactions, const State& state,
                           Kinematics& kinematics ) const
  {
    evaluateKinematics<SequentialExecution>( domain, interactions, state, kinematics );
  }

  template <ExecutionPolicy Execution>
  Real addResidual( const ContactDomain& domain, const InteractionStorage& interactions, const State& state,
                    const Kinematics& kinematics, mfem::Vector& residual ) const
  {
    validateState( state );
    validateResidual( domain, residual );
    const bool use_device = Execution::uses_device;
    const auto view = domain.view( use_device );
    const auto* pairs = interactions.pairs.Read( use_device );
    const auto* pressure = kinematics.pressure.Read( use_device );
    const auto* potential = kinematics.potential_density.Read( use_device );
    const auto* gap = kinematics.gap.Read( use_device );
    auto* output = residual.ReadWrite( use_device );
    const auto parameters = parameters_;
    kinematics.scalar_workspace.UseDevice( use_device );
    kinematics.scalar_workspace = 0.0;
    auto* energy = kinematics.scalar_workspace.ReadWrite( use_device );

#ifdef TRIBOL_USE_ENZYME
    if constexpr ( !ConstraintPolicy::is_nodal && LinearizationPolicy::uses_enzyme ) {
      if ( hasOnlyLinearSegments( domain ) ) {
        forEachInteraction<Execution>( interactions.batches, interactions.coloring,
                                       [=] MFEM_HOST_DEVICE( int pair_index ) {
                                         std::array<Real, max_local_coordinates> coordinates{};
                                         std::array<Real, max_local_coordinates> gradient{};
                                         loadCoordinates( view, pairs[pair_index], coordinates, nullptr );
                                         const EnzymePairContext context{ parameters };
                                         enzymePairGradient( coordinates.data(), &context, gradient.data() );
                                         scatterLocalVector( view, pairs[pair_index], gradient, output );
                                         Real pair_energy{};
                                         enzymePairEnergy( coordinates.data(), &context, &pair_energy );
                                         storeInteractionScalar<Execution>( energy, pair_index, pair_energy );
                                       } );
        finalizeInteractionScalar<Execution>( kinematics.scalar_workspace, interactions.pairs.Size() );
        Execution::synchronize();
        return kinematics.scalar_workspace.HostRead()[0];
      }
    }
#endif
    forEachInteraction<Execution>( interactions.batches, interactions.coloring, [=] MFEM_HOST_DEVICE( int pair_index ) {
      using AD = detail::Gradient<Real, max_local_coordinates>;
      std::array<AD, max_local_coordinates> coordinates{};
      loadCoordinates( view, pairs[pair_index], coordinates, nullptr );
      const auto values = pairValues( view, pairs[pair_index], coordinates, parameters );
      AD pair_energy{};
      if constexpr ( ConstraintPolicy::is_nodal ) {
        const int mortar_element = pairs[pair_index].mortar_element;
        const int mortar_nodes = view.numberOfElementDofs( mortar_element );
        for ( int node = 0; node < mortar_nodes; ++node ) {
          const int dof = view.elementDof( mortar_element, node );
          const Real area_coefficient =
              is_lagrange_multiplier_v<EnforcementPolicy> ? 0.0 : potential[dof] - pressure[dof] * gap[dof];
          pair_energy = pair_energy + AD{ pressure[dof] } * values.weighted_gap[node] +
                        AD{ area_coefficient } * values.area[node];
        }
      } else {
        const auto law = parameters.enforcement;
        for ( int q = 0; q < parameters.quadrature_points; ++q ) {
          pair_energy = pair_energy + values.weights[q] * law.potential( values.gaps[q] );
        }
      }
      scatterGradient( view, pairs[pair_index], pair_energy, output );
      if constexpr ( capabilities.has_energy ) {
        storeInteractionScalar<Execution>( energy, pair_index, detail::scalarValue( pair_energy.value ) );
      }
    } );
    if constexpr ( capabilities.has_energy ) {
      finalizeInteractionScalar<Execution>( kinematics.scalar_workspace, interactions.pairs.Size() );
      Execution::synchronize();
      return kinematics.scalar_workspace.HostRead()[0];
    } else {
      return 0.0;
    }
  }

  Real addResidual( const ContactDomain& domain, const InteractionStorage& interactions, const State& state,
                    const Kinematics& kinematics, mfem::Vector& residual ) const
  {
    return addResidual<SequentialExecution>( domain, interactions, state, kinematics, residual );
  }

  template <ExecutionPolicy Execution>
  void applyCoordinateDerivative( const ContactDomain& domain, const InteractionStorage& interactions,
                                  const State& state, const Kinematics& kinematics, const mfem::Vector& direction,
                                  mfem::Vector& derivative ) const
  {
    validateState( state );
    validateResidual( domain, derivative );
    if ( direction.Size() != domain.coordinateSize() ) {
      throw std::invalid_argument( "EnergyMortar coordinate direction has the wrong size." );
    }
    derivative = 0.0;
    const bool use_device = Execution::uses_device;
    const auto view = domain.view( use_device );
    const auto* pairs = interactions.pairs.Read( use_device );
    const auto* vector = direction.Read( use_device );
    const auto* pressure = kinematics.pressure.Read( use_device );
    const auto* pressure_tangent = kinematics.pressure_tangent.Read( use_device );
    const auto* potential = kinematics.potential_density.Read( use_device );
    const auto* gap = kinematics.gap.Read( use_device );
    auto* output = derivative.ReadWrite( use_device );
    const auto parameters = parameters_;

#ifdef TRIBOL_USE_ENZYME
    if constexpr ( !ConstraintPolicy::is_nodal && LinearizationPolicy::uses_enzyme ) {
      if ( hasOnlyLinearSegments( domain ) ) {
        forEachInteraction<Execution>(
            interactions.batches, interactions.coloring, [=] MFEM_HOST_DEVICE( int pair_index ) {
              std::array<Real, max_local_coordinates> coordinates{};
              std::array<Real, max_local_coordinates> local_direction{};
              std::array<Real, max_local_coordinates> gradient{};
              std::array<Real, max_local_coordinates> derivative_values{};
              loadCoordinates( view, pairs[pair_index], coordinates, nullptr );
              loadDirection( view, pairs[pair_index], vector, local_direction );
              const EnzymePairContext context{ parameters };
              __enzyme_fwddiff<void>( (void*)enzymePairGradient, TRIBOL_FUTURE_ENZYME_DUP, coordinates.data(),
                                      local_direction.data(), TRIBOL_FUTURE_ENZYME_CONST, (const void*)&context,
                                      TRIBOL_FUTURE_ENZYME_DUP, gradient.data(), derivative_values.data() );
              scatterLocalVector( view, pairs[pair_index], derivative_values, output );
            } );
        return;
      }
    }
#endif

    kinematics.directional_weighted_gap.UseDevice( use_device );
    kinematics.directional_area.UseDevice( use_device );
    kinematics.directional_weighted_gap = 0.0;
    kinematics.directional_area = 0.0;
    auto* dweighted_gap = kinematics.directional_weighted_gap.ReadWrite( use_device );
    auto* darea = kinematics.directional_area.ReadWrite( use_device );
    const auto* area_values = kinematics.tributary_area.Read( use_device );
    const auto* weighted_values = kinematics.weighted_gap.Read( use_device );

    if constexpr ( ConstraintPolicy::is_nodal && !is_lagrange_multiplier_v<EnforcementPolicy> ) {
      forEachInteraction<Execution>(
          interactions.batches, interactions.coloring, [=] MFEM_HOST_DEVICE( int pair_index ) {
            using Directional = detail::Gradient<Real, 1>;
            std::array<Directional, max_local_coordinates> coordinates{};
            loadCoordinates( view, pairs[pair_index], coordinates, vector );
            const auto values = pairValues( view, pairs[pair_index], coordinates, parameters );
            const int mortar_element = pairs[pair_index].mortar_element;
            const int mortar_nodes = view.numberOfElementDofs( mortar_element );
            for ( int node = 0; node < mortar_nodes; ++node ) {
              const int dof = view.elementDof( mortar_element, node );
              AtomicAdd( dweighted_gap[dof], values.weighted_gap[node].derivative[0] );
              AtomicAdd( darea[dof], values.area[node].derivative[0] );
            }
          } );
    }

    forEachInteraction<Execution>( interactions.batches, interactions.coloring, [=] MFEM_HOST_DEVICE( int pair_index ) {
      using Directional = detail::Gradient<Real, 1>;
      using AD = detail::Gradient<Directional, max_local_coordinates>;
      std::array<AD, max_local_coordinates> coordinates{};
      loadNestedCoordinates( view, pairs[pair_index], coordinates, vector );
      const auto values = pairValues( view, pairs[pair_index], coordinates, parameters );
      AD pair_energy{};
      if constexpr ( ConstraintPolicy::is_nodal ) {
        const int mortar_element = pairs[pair_index].mortar_element;
        const int mortar_nodes = view.numberOfElementDofs( mortar_element );
        for ( int node = 0; node < mortar_nodes; ++node ) {
          const int dof = view.elementDof( mortar_element, node );
          Real pressure_dot = 0.0;
          Real area_coefficient = 0.0;
          Real area_coefficient_dot = 0.0;
          if constexpr ( !is_lagrange_multiplier_v<EnforcementPolicy> ) {
            const Real area_value = area_values[dof];
            const Real weighted_value = weighted_values[dof];
            const Real gap_dot = area_value > 0.0 ? ( dweighted_gap[dof] * area_value - weighted_value * darea[dof] ) /
                                                        ( area_value * area_value )
                                                  : 0.0;
            pressure_dot = pressure_tangent[dof] * gap_dot;
            area_coefficient = potential[dof] - pressure[dof] * gap[dof];
            area_coefficient_dot = -gap[dof] * pressure_dot;
          }
          const AD pressure_ad{ Directional::variable( pressure[dof], 0 ) };
          AD corrected_pressure = pressure_ad;
          corrected_pressure.value.derivative[0] = pressure_dot;
          AD area_coefficient_ad{ Directional{ area_coefficient } };
          area_coefficient_ad.value.derivative[0] = area_coefficient_dot;
          pair_energy =
              pair_energy + corrected_pressure * values.weighted_gap[node] + area_coefficient_ad * values.area[node];
        }
      } else {
        const auto law = parameters.enforcement;
        for ( int q = 0; q < parameters.quadrature_points; ++q ) {
          pair_energy = pair_energy + values.weights[q] * law.potential( values.gaps[q] );
        }
      }
      scatterHessianVector( view, pairs[pair_index], pair_energy, output );
    } );
  }

  void applyCoordinateDerivative( const ContactDomain& domain, const InteractionStorage& interactions,
                                  const State& state, const Kinematics& kinematics, const mfem::Vector& direction,
                                  mfem::Vector& derivative ) const
  {
    applyCoordinateDerivative<SequentialExecution>( domain, interactions, state, kinematics, direction, derivative );
  }

  template <ExecutionPolicy Execution>
  void applyGapDerivative( const ContactDomain& domain, const InteractionStorage& interactions,
                           const Kinematics& kinematics, const mfem::Vector& direction,
                           mfem::Vector& gap_direction ) const
  {
    if ( direction.Size() != domain.coordinateSize() || gap_direction.Size() != domain.numberOfDofs() ) {
      throw std::invalid_argument( "EnergyMortar gap derivative received an incorrectly sized vector." );
    }
    const bool use_device = Execution::uses_device;
    kinematics.directional_area.UseDevice( use_device );
    gap_direction.UseDevice( use_device );
    gap_direction = 0.0;
    kinematics.directional_area = 0.0;
    auto* weighted_dot = gap_direction.ReadWrite( use_device );
    auto* area_dot = kinematics.directional_area.ReadWrite( use_device );
    const auto view = domain.view( use_device );
    const auto* pairs = interactions.pairs.Read( use_device );
    const auto* vector = direction.Read( use_device );
    const auto parameters = parameters_;
    forEachInteraction<Execution>( interactions.batches, interactions.coloring, [=] MFEM_HOST_DEVICE( int pair_index ) {
      using Directional = detail::Gradient<Real, 1>;
      std::array<Directional, max_local_coordinates> coordinates{};
      loadCoordinates( view, pairs[pair_index], coordinates, vector );
      const auto values = pairValues( view, pairs[pair_index], coordinates, parameters );
      const int mortar_element = pairs[pair_index].mortar_element;
      const int mortar_nodes = view.numberOfElementDofs( mortar_element );
      for ( int node = 0; node < mortar_nodes; ++node ) {
        const int dof = view.elementDof( mortar_element, node );
        AtomicAdd( weighted_dot[dof], values.weighted_gap[node].derivative[0] );
        AtomicAdd( area_dot[dof], values.area[node].derivative[0] );
      }
    } );
  }

  void applyGapDerivative( const ContactDomain& domain, const InteractionStorage& interactions,
                           const Kinematics& kinematics, const mfem::Vector& direction,
                           mfem::Vector& gap_direction ) const
  {
    applyGapDerivative<SequentialExecution>( domain, interactions, kinematics, direction, gap_direction );
  }

  template <ExecutionPolicy Execution>
  void applyDualDerivative( const ContactDomain& domain, const InteractionStorage& interactions,
                            const mfem::Vector& dual_direction, mfem::Vector& force_direction ) const
  {
    static_assert( is_lagrange_multiplier_v<EnforcementPolicy>,
                   "applyDualDerivative is available only for Lagrange-multiplier EnergyMortar." );
    if ( dual_direction.Size() != domain.numberOfDofs() ) {
      throw std::invalid_argument( "EnergyMortar dual direction has the wrong size." );
    }
    validateResidual( domain, force_direction );
    force_direction = 0.0;
    const bool use_device = Execution::uses_device;
    const auto view = domain.view( use_device );
    const auto* pairs = interactions.pairs.Read( use_device );
    const auto* dual = dual_direction.Read( use_device );
    auto* output = force_direction.ReadWrite( use_device );
    const auto parameters = parameters_;
    forEachInteraction<Execution>( interactions.batches, interactions.coloring, [=] MFEM_HOST_DEVICE( int pair_index ) {
      using AD = detail::Gradient<Real, max_local_coordinates>;
      std::array<AD, max_local_coordinates> coordinates{};
      loadCoordinates( view, pairs[pair_index], coordinates, nullptr );
      const auto values = pairValues( view, pairs[pair_index], coordinates, parameters );
      AD pair_energy{};
      const int mortar_element = pairs[pair_index].mortar_element;
      const int mortar_nodes = view.numberOfElementDofs( mortar_element );
      for ( int node = 0; node < mortar_nodes; ++node ) {
        pair_energy = pair_energy + AD{ dual[view.elementDof( mortar_element, node )] } * values.weighted_gap[node];
      }
      scatterGradient( view, pairs[pair_index], pair_energy, output );
    } );
  }

  void applyDualDerivative( const ContactDomain& domain, const InteractionStorage& interactions,
                            const mfem::Vector& dual_direction, mfem::Vector& force_direction ) const
  {
    applyDualDerivative<SequentialExecution>( domain, interactions, dual_direction, force_direction );
  }

  const Parameters& parameters() const { return parameters_; }

 private:
  static bool hasOnlyLinearSegments( const ContactDomain& domain )
  {
    for ( int element = 0; element < domain.elementOrders().Size(); ++element ) {
      if ( domain.elementOrders()[element] != 1 ) {
        return false;
      }
    }
    return true;
  }

#ifdef TRIBOL_USE_ENZYME
  struct EnzymePairContext {
    Parameters parameters;
  };

  TRIBOL_FUTURE_HOST_DEVICE static void enzymePairEnergy( const Real* coordinates, const void* context_pointer,
                                                          Real* energy )
  {
    const auto& parameters = static_cast<const EnzymePairContext*>( context_pointer )->parameters;
    const int nonmortar_offset = 2 * max_nodes_per_segment;
    const Real mortar_x = coordinates[2] - coordinates[0];
    const Real mortar_y = coordinates[3] - coordinates[1];
    const Real nonmortar_x = coordinates[nonmortar_offset + 2] - coordinates[nonmortar_offset];
    const Real nonmortar_y = coordinates[nonmortar_offset + 3] - coordinates[nonmortar_offset + 1];
    const Real nonmortar_length = std::sqrt( nonmortar_x * nonmortar_x + nonmortar_y * nonmortar_y );
    const Real normal_x = nonmortar_y / nonmortar_length;
    const Real normal_y = -nonmortar_x / nonmortar_length;
    const Real mortar_denominator = mortar_x * mortar_x + mortar_y * mortar_y;
    Real lower = ( ( coordinates[nonmortar_offset] - coordinates[0] ) * mortar_x +
                   ( coordinates[nonmortar_offset + 1] - coordinates[1] ) * mortar_y ) /
                 mortar_denominator;
    Real upper = ( ( coordinates[nonmortar_offset + 2] - coordinates[0] ) * mortar_x +
                   ( coordinates[nonmortar_offset + 3] - coordinates[1] ) * mortar_y ) /
                 mortar_denominator;
    if ( upper < lower ) {
      const Real temporary = lower;
      lower = upper;
      upper = temporary;
    }
    lower = lower < 0.0 ? 0.0 : ( lower > 1.0 ? 1.0 : lower );
    upper = upper < 0.0 ? 0.0 : ( upper > 1.0 ? 1.0 : upper );
    const Real half_interval = 0.5 * ( upper - lower );
    const Real midpoint = 0.5 * ( upper + lower );
    Real result{};
    for ( int quadrature_point = 0; quadrature_point < parameters.quadrature_points; ++quadrature_point ) {
      Real point{};
      Real weight{};
      gaussPoint( parameters.quadrature_points, quadrature_point, point, weight );
      const Real xi = midpoint + half_interval * point;
      const Real mortar_point_x = ( 1.0 - xi ) * coordinates[0] + xi * coordinates[2];
      const Real mortar_point_y = ( 1.0 - xi ) * coordinates[1] + xi * coordinates[3];
      const Real eta = ( ( mortar_point_x - coordinates[nonmortar_offset] ) * nonmortar_x +
                         ( mortar_point_y - coordinates[nonmortar_offset + 1] ) * nonmortar_y ) /
                       ( nonmortar_x * nonmortar_x + nonmortar_y * nonmortar_y );
      const Real nonmortar_point_x = coordinates[nonmortar_offset] + eta * nonmortar_x;
      const Real nonmortar_point_y = coordinates[nonmortar_offset + 1] + eta * nonmortar_y;
      Real gap = ( mortar_point_x - nonmortar_point_x ) * normal_x + ( mortar_point_y - nonmortar_point_y ) * normal_y;
      if ( parameters.endpoint_smoothing > 0.0 ) {
        const Real distance = 2.0 * std::min( xi, 1.0 - xi );
        gap *= std::min( 1.0, distance / parameters.endpoint_smoothing );
      }
      const Real mortar_length = std::sqrt( mortar_x * mortar_x + mortar_y * mortar_y );
      const Real physical_weight = weight * half_interval * mortar_length;
      result += physical_weight * parameters.enforcement.potential( gap );
    }
    *energy = result;
  }

  TRIBOL_FUTURE_HOST_DEVICE static void enzymePairGradient( const Real* coordinates, const void* context,
                                                            Real* gradient )
  {
    Real coordinate_shadow[max_local_coordinates]{};
    Real energy{};
    Real energy_shadow{ 1.0 };
    __enzyme_autodiff<void>( (void*)enzymePairEnergy, TRIBOL_FUTURE_ENZYME_DUP, coordinates, coordinate_shadow,
                             TRIBOL_FUTURE_ENZYME_CONST, context, TRIBOL_FUTURE_ENZYME_DUP, &energy, &energy_shadow );
    for ( int coordinate = 0; coordinate < max_local_coordinates; ++coordinate ) {
      gradient[coordinate] = coordinate_shadow[coordinate];
    }
  }
#endif

  template <typename T>
  struct PairValues {
    std::array<T, max_nodes_per_segment> weighted_gap{};
    std::array<T, max_nodes_per_segment> area{};
    std::array<T, 3> gaps{};
    std::array<T, 3> weights{};
    std::array<detail::Vector<T, 2>, 3> mortar_points{};
    std::array<detail::Vector<T, 2>, 3> nonmortar_points{};
  };

  template <typename T>
  struct CurveValues {
    detail::Vector<T, 2> point{};
    detail::Vector<T, 2> tangent{};
    detail::Vector<T, 2> second_derivative{};
    std::array<T, max_nodes_per_segment> shape{};
  };

  TRIBOL_FUTURE_HOST_DEVICE static void referenceRange( const ContactDomainView& domain, int element, Real& lower,
                                                        Real& upper )
  {
    lower = domain.referenceCoordinate( element, 0 );
    upper = lower;
    for ( int node = 1; node < domain.numberOfElementDofs( element ); ++node ) {
      const Real coordinate = domain.referenceCoordinate( element, node );
      lower = coordinate < lower ? coordinate : lower;
      upper = coordinate > upper ? coordinate : upper;
    }
  }

  template <typename T>
  TRIBOL_FUTURE_HOST_DEVICE static CurveValues<T> evaluateCurve(
      const ContactDomainView& domain, int element, int side, const std::array<T, max_local_coordinates>& coordinates,
      const T& parameter )
  {
    using detail::operator+;
    using detail::operator*;
    CurveValues<T> result;
    const int number_of_nodes = domain.numberOfElementDofs( element );
    for ( int node = 0; node < number_of_nodes; ++node ) {
      const Real node_coordinate = domain.referenceCoordinate( element, node );
      T shape{ 1.0 };
      for ( int other = 0; other < number_of_nodes; ++other ) {
        if ( other != node ) {
          const Real other_coordinate = domain.referenceCoordinate( element, other );
          shape = shape * ( parameter - T{ other_coordinate } ) / T{ node_coordinate - other_coordinate };
        }
      }
      result.shape[node] = shape;

      T derivative{};
      for ( int omitted = 0; omitted < number_of_nodes; ++omitted ) {
        if ( omitted == node ) {
          continue;
        }
        const Real omitted_coordinate = domain.referenceCoordinate( element, omitted );
        T term{ 1.0 / ( node_coordinate - omitted_coordinate ) };
        for ( int other = 0; other < number_of_nodes; ++other ) {
          if ( other != node && other != omitted ) {
            const Real other_coordinate = domain.referenceCoordinate( element, other );
            term = term * ( parameter - T{ other_coordinate } ) / T{ node_coordinate - other_coordinate };
          }
        }
        derivative = derivative + term;
      }

      T second_derivative{};
      for ( int first_omitted = 0; first_omitted < number_of_nodes; ++first_omitted ) {
        if ( first_omitted == node ) {
          continue;
        }
        for ( int second_omitted = 0; second_omitted < number_of_nodes; ++second_omitted ) {
          if ( second_omitted == node || second_omitted == first_omitted ) {
            continue;
          }
          T term{ 1.0 / ( ( node_coordinate - domain.referenceCoordinate( element, first_omitted ) ) *
                          ( node_coordinate - domain.referenceCoordinate( element, second_omitted ) ) ) };
          for ( int other = 0; other < number_of_nodes; ++other ) {
            if ( other != node && other != first_omitted && other != second_omitted ) {
              const Real other_coordinate = domain.referenceCoordinate( element, other );
              term = term * ( parameter - T{ other_coordinate } ) / T{ node_coordinate - other_coordinate };
            }
          }
          second_derivative = second_derivative + term;
        }
      }

      detail::Vector<T, 2> point{};
      for ( int component = 0; component < 2; ++component ) {
        point[component] = coordinates[( side * max_nodes_per_segment + node ) * 2 + component];
      }
      result.point = result.point + shape * point;
      result.tangent = result.tangent + derivative * point;
      result.second_derivative = result.second_derivative + second_derivative * point;
    }
    return result;
  }

  template <typename T>
  TRIBOL_FUTURE_HOST_DEVICE static T projectPointToCurve( const ContactDomainView& domain, int element, int side,
                                                          const std::array<T, max_local_coordinates>& coordinates,
                                                          const detail::Vector<T, 2>& target,
                                                          const Parameters& parameters )
  {
    using detail::operator-;
    Real lower{};
    Real upper{};
    referenceRange( domain, element, lower, upper );
    const auto first = evaluateCurve( domain, element, side, coordinates, T{ lower } );
    const auto last = evaluateCurve( domain, element, side, coordinates, T{ upper } );
    const auto chord = last.point - first.point;
    const T chord_norm = detail::dot( chord, chord );
    T parameter{ 0.5 * ( lower + upper ) };
    if ( std::abs( detail::scalarValue( chord_norm ) ) > 1.0e-28 ) {
      parameter = T{ lower } + T{ upper - lower } * detail::dot( target - first.point, chord ) / chord_norm;
      parameter = detail::clamp( parameter, lower, upper );
    }
    for ( int iteration = 0; iteration < parameters.projection_iterations; ++iteration ) {
      const auto curve = evaluateCurve( domain, element, side, coordinates, parameter );
      const auto difference = curve.point - target;
      const T numerator = detail::dot( difference, curve.tangent );
      const T denominator =
          detail::dot( curve.tangent, curve.tangent ) + detail::dot( difference, curve.second_derivative );
      if ( std::abs( detail::scalarValue( denominator ) ) > 1.0e-28 ) {
        parameter = detail::clamp( parameter - numerator / denominator, lower, upper );
      }
    }
    return parameter;
  }

  template <typename T>
  TRIBOL_FUTURE_HOST_DEVICE static PairValues<T> pairValues( const ContactDomainView& domain, ElementPair pair,
                                                             const std::array<T, max_local_coordinates>& coordinates,
                                                             const Parameters& parameters )
  {
    using detail::operator-;
    const int mortar = pair.mortar_element;
    const int nonmortar = pair.nonmortar_element;
    Real mortar_lower{};
    Real mortar_upper{};
    Real nonmortar_lower{};
    Real nonmortar_upper{};
    referenceRange( domain, mortar, mortar_lower, mortar_upper );
    referenceRange( domain, nonmortar, nonmortar_lower, nonmortar_upper );
    const auto nonmortar_first = evaluateCurve( domain, nonmortar, 1, coordinates, T{ nonmortar_lower } );
    const auto nonmortar_last = evaluateCurve( domain, nonmortar, 1, coordinates, T{ nonmortar_upper } );
    T lower = projectPointToCurve( domain, mortar, 0, coordinates, nonmortar_first.point, parameters );
    T upper = projectPointToCurve( domain, mortar, 0, coordinates, nonmortar_last.point, parameters );
    if ( detail::scalarValue( upper ) < detail::scalarValue( lower ) ) {
      const T temporary = lower;
      lower = upper;
      upper = temporary;
    }
    lower = detail::clamp( lower, mortar_lower, mortar_upper );
    upper = detail::clamp( upper, mortar_lower, mortar_upper );

    PairValues<T> result;
    const T half_interval = T{ 0.5 } * ( upper - lower );
    const T midpoint = T{ 0.5 } * ( upper + lower );
    for ( int quadrature_point = 0; quadrature_point < parameters.quadrature_points; ++quadrature_point ) {
      Real point{};
      Real weight{};
      gaussPoint( parameters.quadrature_points, quadrature_point, point, weight );
      const T xi = midpoint + half_interval * T{ point };
      const auto mortar_curve = evaluateCurve( domain, mortar, 0, coordinates, xi );
      const T eta = projectPointToCurve( domain, nonmortar, 1, coordinates, mortar_curve.point, parameters );
      const auto nonmortar_curve = evaluateCurve( domain, nonmortar, 1, coordinates, eta );
      const T nonmortar_length = detail::norm( nonmortar_curve.tangent );
      if ( std::abs( detail::scalarValue( nonmortar_length ) ) <= 1.0e-28 ) {
        continue;
      }
      const detail::Vector<T, 2> normal{ nonmortar_curve.tangent[1] / nonmortar_length,
                                         -nonmortar_curve.tangent[0] / nonmortar_length };
      T gap = detail::dot( mortar_curve.point - nonmortar_curve.point, normal );
      if ( parameters.endpoint_smoothing > 0.0 ) {
        const Real xi_value = detail::scalarValue( xi );
        const Real normalized = ( xi_value - mortar_lower ) / ( mortar_upper - mortar_lower );
        const Real distance = 2.0 * std::min( normalized, 1.0 - normalized );
        const Real smoothing = std::min( 1.0, distance / parameters.endpoint_smoothing );
        gap = T{ smoothing } * gap;
      }
      const T physical_weight = T{ weight } * half_interval * detail::norm( mortar_curve.tangent );
      result.gaps[quadrature_point] = gap;
      result.weights[quadrature_point] = physical_weight;
      result.mortar_points[quadrature_point] = mortar_curve.point;
      result.nonmortar_points[quadrature_point] = nonmortar_curve.point;
      for ( int node = 0; node < domain.numberOfElementDofs( mortar ); ++node ) {
        result.weighted_gap[node] = result.weighted_gap[node] + physical_weight * mortar_curve.shape[node] * gap;
        result.area[node] = result.area[node] + physical_weight * mortar_curve.shape[node];
      }
    }
    return result;
  }

  template <typename T>
  TRIBOL_FUTURE_HOST_DEVICE static PairValues<T> pairValues( const ContactDomainView& domain, ElementPair pair,
                                                             const Parameters& parameters )
  {
    std::array<T, max_local_coordinates> coordinates{};
    loadCoordinates( domain, pair, coordinates, nullptr );
    return pairValues( domain, pair, coordinates, parameters );
  }

  TRIBOL_FUTURE_HOST_DEVICE static void gaussPoint( int number_of_points, int point_index, Real& point, Real& weight )
  {
    if ( number_of_points == 1 ) {
      point = 0.0;
      weight = 2.0;
    } else if ( number_of_points == 2 ) {
      constexpr Real coordinate = 0.57735026918962576451;
      point = point_index == 0 ? -coordinate : coordinate;
      weight = 1.0;
    } else {
      constexpr Real coordinate = 0.77459666924148337704;
      point = point_index == 0 ? -coordinate : ( point_index == 1 ? 0.0 : coordinate );
      weight = point_index == 1 ? 8.0 / 9.0 : 5.0 / 9.0;
    }
  }

  Real overlapLength( const ContactDomainView& domain, ElementPair pair ) const
  {
    auto overlap_parameters = parameters_;
    overlap_parameters.quadrature_points = 1;
    const auto values = pairValues<Real>( domain, pair, overlap_parameters );
    Real length{};
    for ( int q = 0; q < overlap_parameters.quadrature_points; ++q ) {
      length += values.weights[q];
    }
    return std::abs( length );
  }

  template <typename T>
  TRIBOL_FUTURE_HOST_DEVICE static void loadCoordinates( const ContactDomainView& domain, ElementPair pair,
                                                         std::array<T, max_local_coordinates>& coordinates,
                                                         const Real* direction )
  {
    const int elements[2] = { pair.mortar_element, pair.nonmortar_element };
    for ( int side = 0; side < 2; ++side ) {
      const int number_of_nodes = domain.numberOfElementDofs( elements[side] );
      for ( int node = 0; node < number_of_nodes; ++node ) {
        const int dof = domain.elementDof( elements[side], node );
        for ( int component = 0; component < 2; ++component ) {
          const int local = ( side * max_nodes_per_segment + node ) * 2 + component;
          const int global = component * domain.number_of_dofs + dof;
          if constexpr ( std::is_same_v<T, Real> ) {
            coordinates[local] = domain.coordinates[global];
          } else {
            if ( direction != nullptr ) {
              coordinates[local] = T{ domain.coordinates[global] };
              coordinates[local].derivative[0] = direction[global];
            } else {
              coordinates[local] = T::variable( domain.coordinates[global], local );
            }
          }
        }
      }
    }
  }

  template <typename T>
  TRIBOL_FUTURE_HOST_DEVICE static void loadNestedCoordinates( const ContactDomainView& domain, ElementPair pair,
                                                               std::array<T, max_local_coordinates>& coordinates,
                                                               const Real* direction )
  {
    using Directional = decltype( T{}.value );
    const int elements[2] = { pair.mortar_element, pair.nonmortar_element };
    for ( int side = 0; side < 2; ++side ) {
      const int number_of_nodes = domain.numberOfElementDofs( elements[side] );
      for ( int node = 0; node < number_of_nodes; ++node ) {
        const int dof = domain.elementDof( elements[side], node );
        for ( int component = 0; component < 2; ++component ) {
          const int local = ( side * max_nodes_per_segment + node ) * 2 + component;
          const int global = component * domain.number_of_dofs + dof;
          Directional value( domain.coordinates[global] );
          value.derivative[0] = direction[global];
          coordinates[local] = T::variable( value, local );
        }
      }
    }
  }

  TRIBOL_FUTURE_HOST_DEVICE static void loadDirection( const ContactDomainView& domain, ElementPair pair,
                                                       const Real* direction,
                                                       std::array<Real, max_local_coordinates>& local_direction )
  {
    const int elements[2] = { pair.mortar_element, pair.nonmortar_element };
    for ( int side = 0; side < 2; ++side ) {
      const int number_of_nodes = domain.numberOfElementDofs( elements[side] );
      for ( int node = 0; node < number_of_nodes; ++node ) {
        const int dof = domain.elementDof( elements[side], node );
        for ( int component = 0; component < 2; ++component ) {
          local_direction[( side * max_nodes_per_segment + node ) * 2 + component] =
              direction[component * domain.number_of_dofs + dof];
        }
      }
    }
  }

  template <typename T>
  TRIBOL_FUTURE_HOST_DEVICE static void scatterLocalVector( const ContactDomainView& domain, ElementPair pair,
                                                            const std::array<T, max_local_coordinates>& local_values,
                                                            Real* output )
  {
    const int elements[2] = { pair.mortar_element, pair.nonmortar_element };
    for ( int side = 0; side < 2; ++side ) {
      const int number_of_nodes = domain.numberOfElementDofs( elements[side] );
      for ( int node = 0; node < number_of_nodes; ++node ) {
        const int dof = domain.elementDof( elements[side], node );
        for ( int component = 0; component < 2; ++component ) {
          const int local = ( side * max_nodes_per_segment + node ) * 2 + component;
          AtomicAdd( output[component * domain.number_of_dofs + dof], detail::scalarValue( local_values[local] ) );
        }
      }
    }
  }

  template <typename AD>
  TRIBOL_FUTURE_HOST_DEVICE static void scatterGradient( const ContactDomainView& domain, ElementPair pair,
                                                         const AD& energy, Real* residual )
  {
    const int elements[2] = { pair.mortar_element, pair.nonmortar_element };
    for ( int side = 0; side < 2; ++side ) {
      const int number_of_nodes = domain.numberOfElementDofs( elements[side] );
      for ( int node = 0; node < number_of_nodes; ++node ) {
        const int dof = domain.elementDof( elements[side], node );
        for ( int component = 0; component < 2; ++component ) {
          const int local = ( side * max_nodes_per_segment + node ) * 2 + component;
          AtomicAdd( residual[component * domain.number_of_dofs + dof],
                     detail::scalarValue( energy.derivative[local] ) );
        }
      }
    }
  }

  template <typename AD>
  TRIBOL_FUTURE_HOST_DEVICE static void scatterHessianVector( const ContactDomainView& domain, ElementPair pair,
                                                              const AD& energy, Real* derivative )
  {
    const int elements[2] = { pair.mortar_element, pair.nonmortar_element };
    for ( int side = 0; side < 2; ++side ) {
      const int number_of_nodes = domain.numberOfElementDofs( elements[side] );
      for ( int node = 0; node < number_of_nodes; ++node ) {
        const int dof = domain.elementDof( elements[side], node );
        for ( int component = 0; component < 2; ++component ) {
          const int local = ( side * max_nodes_per_segment + node ) * 2 + component;
          AtomicAdd( derivative[component * domain.number_of_dofs + dof], energy.derivative[local].derivative[0] );
        }
      }
    }
  }

  static void validateResidual( const ContactDomain& domain, const mfem::Vector& residual )
  {
    if ( residual.Size() != domain.coordinateSize() ) {
      throw std::invalid_argument( "EnergyMortar residual has the wrong size." );
    }
  }

  void validateState( const State& state ) const
  {
    if constexpr ( is_lagrange_multiplier_v<EnforcementPolicy> ) {
      if ( state.dual == nullptr ) {
        throw std::invalid_argument( "Nodal Lagrange-multiplier EnergyMortar requires a dual vector." );
      }
    }
    if constexpr ( is_external_pressure_v<EnforcementPolicy> ) {
      if ( state.external_pressure == nullptr ) {
        throw std::invalid_argument( "External-pressure EnergyMortar requires staged pressure data." );
      }
      const int expected_size = state.external_pressure->pressure.Size();
      if ( state.external_pressure->potential_density.Size() != expected_size ||
           state.external_pressure->pressure_tangent.Size() != expected_size ) {
        throw std::invalid_argument( "External pressure, tangent, and potential arrays must have equal sizes." );
      }
    }
  }

  static void setMemoryMode( Kinematics& kinematics, const InteractionStorage& interactions, bool use_device )
  {
    kinematics.gap.UseDevice( use_device );
    kinematics.weighted_gap.UseDevice( use_device );
    kinematics.tributary_area.UseDevice( use_device );
    kinematics.pressure.UseDevice( use_device );
    kinematics.pressure_tangent.UseDevice( use_device );
    kinematics.potential_density.UseDevice( use_device );
    kinematics.pair_offsets.GetMemory().UseDevice( use_device );
    kinematics.quadrature_points.UseDevice( use_device );
    kinematics.quadrature_weights.UseDevice( use_device );
    kinematics.quadrature_gaps.UseDevice( use_device );
    kinematics.directional_weighted_gap.UseDevice( use_device );
    kinematics.directional_area.UseDevice( use_device );
    kinematics.scalar_workspace.UseDevice( use_device );
    const_cast<mfem::Array<ElementPair>&>( interactions.pairs ).GetMemory().UseDevice( use_device );
    const_cast<InteractionColoring&>( interactions.coloring ).offsets.GetMemory().UseDevice( use_device );
    const_cast<InteractionColoring&>( interactions.coloring ).pair_indices.GetMemory().UseDevice( use_device );
  }

  Parameters parameters_;
};

}  // namespace tribol::future

#endif
