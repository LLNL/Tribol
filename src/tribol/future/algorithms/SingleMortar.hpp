#ifndef TRIBOL_FUTURE_ALGORITHMS_SINGLEMORTAR_HPP_
#define TRIBOL_FUTURE_ALGORITHMS_SINGLEMORTAR_HPP_

#include "tribol/future/Concepts.hpp"
#include "tribol/future/ContactDomain.hpp"
#include "tribol/future/Execution.hpp"
#include "tribol/future/InteractionBatches.hpp"
#include "tribol/future/Policies.hpp"
#include "tribol/future/Results.hpp"
#include "tribol/future/detail/AutomaticDifferentiation.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace tribol::future {

template <MortarBasisPolicy BasisPolicy = DualMortarBasis, typename EnforcementPolicy = NodalLagrangeMultiplier,
          LinearizationPolicy LinearizationPolicy = AnalyticLinearization>
class SingleMortar {
 public:
  static constexpr int max_nodes_per_element = 4;
  static constexpr int max_local_coordinates = 2 * max_nodes_per_element * 3;
  static constexpr int maximum_element_order = 1;

  struct Parameters {
    BasisPolicy basis{};
    int quadrature_points{ 3 };
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
    mutable mfem::Vector directional_area;
  };

  struct InteractionStorage {
    mfem::Array<ElementPair> pairs;
    mfem::Array<InteractionBatch> batches;
    InteractionColoring coloring;
  };

  struct Linearization {
    using policy_type = LinearizationPolicy;
  };

  static constexpr bool accepts_external_pressure = false;
  static constexpr bool has_nodal_kinematics = true;
  static constexpr MethodCapabilities capabilities{
      .needs_dual_space = true,
      .needs_velocity = false,
      .needs_reference_field = false,
      .needs_material_fields = false,
      .has_history = false,
      .has_energy = false,
      .supports_self_contact = false,
      .supports_assembled_jacobian = true,
  };

  explicit SingleMortar( Parameters parameters = {} ) : parameters_( parameters )
  {
    static_assert( is_lagrange_multiplier_v<EnforcementPolicy>,
                   "Classical SingleMortar currently requires nodal Lagrange-multiplier enforcement." );
    if ( parameters_.quadrature_points != 1 && parameters_.quadrature_points != 3 ) {
      throw std::invalid_argument( "SingleMortar supports one- or three-point triangle quadrature." );
    }
  }

  void buildInteractions( const ContactDomain& domain, const CandidatePairs& candidates,
                          InteractionStorage& interactions ) const
  {
    const auto view = domain.view( false );
    if ( view.dimension != 3 ) {
      throw std::invalid_argument( "SingleMortar supports three-dimensional surface contact." );
    }
    std::vector<ElementPair> accepted;
    accepted.reserve( candidates.pairs.Size() );
    for ( int pair_index = 0; pair_index < candidates.pairs.Size(); ++pair_index ) {
      const auto pair = candidates.pairs[pair_index];
      if ( view.orders[pair.mortar_element] != 1 || view.orders[pair.nonmortar_element] != 1 ) {
        throw std::invalid_argument( "SingleMortar requires linear surface elements; use LowOrderRefinedSurface." );
      }
      if ( isFace( view.topology( pair.mortar_element ) ) && isFace( view.topology( pair.nonmortar_element ) ) &&
           detail::scalarValue( pairValues<Real>( view, pair, nullptr, parameters_ ).measure ) >
               parameters_.minimum_overlap ) {
        accepted.push_back( pair );
      }
    }
    buildInteractionBatches( domain, std::move( accepted ), interactions.pairs, interactions.batches );
    buildInteractionColoring( domain, interactions.pairs, interactions.coloring );
  }

  void prepare( const ContactDomain& domain, const InteractionStorage& interactions, Kinematics& kinematics ) const
  {
    const int dofs = domain.numberOfDofs();
    const int pairs = interactions.pairs.Size();
    kinematics.gap.SetSize( dofs );
    kinematics.weighted_gap.SetSize( dofs );
    kinematics.tributary_area.SetSize( dofs );
    kinematics.pressure.SetSize( dofs );
    kinematics.pressure_tangent.SetSize( dofs );
    kinematics.potential_density.SetSize( dofs );
    kinematics.directional_area.SetSize( dofs );
    kinematics.pair_offsets.SetSize( pairs + 1 );
    kinematics.quadrature_points.SetSize( pairs * 6 );
    kinematics.quadrature_weights.SetSize( pairs );
    kinematics.quadrature_gaps.SetSize( pairs );
    for ( int pair = 0; pair <= pairs; ++pair ) {
      kinematics.pair_offsets[pair] = pair;
    }
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
    const auto view = domain.view( use_device );
    const auto* pairs = interactions.pairs.Read( use_device );
    auto* weighted_gap = kinematics.weighted_gap.ReadWrite( use_device );
    auto* area = kinematics.tributary_area.ReadWrite( use_device );
    auto* diagnostic_gap = kinematics.quadrature_gaps.Write( use_device );
    auto* diagnostic_weight = kinematics.quadrature_weights.Write( use_device );
    auto* diagnostic_points = kinematics.quadrature_points.Write( use_device );
    const Parameters parameters = parameters_;
    forEachInteraction<Execution>( interactions.batches, interactions.coloring, [=] MFEM_HOST_DEVICE( int pair_index ) {
      const auto geometry = pairValues<Real>( view, pairs[pair_index], nullptr, parameters );
      const int mortar = pairs[pair_index].mortar_element;
      const int nodes = std::min( view.numberOfElementDofs( mortar ), max_nodes_per_element );
      for ( int node = 0; node < nodes; ++node ) {
        const int dof = view.elementDof( mortar, node );
        AtomicAdd( weighted_gap[dof], geometry.weighted_gap[node] );
        AtomicAdd( area[dof], geometry.area[node] );
      }
      diagnostic_gap[pair_index] = geometry.gap;
      diagnostic_weight[pair_index] = geometry.measure;
      for ( int component = 0; component < 3; ++component ) {
        diagnostic_points[6 * pair_index + component] = geometry.mortar_centroid[component];
        diagnostic_points[6 * pair_index + 3 + component] = geometry.nonmortar_centroid[component];
      }
    } );
    auto* gap = kinematics.gap.ReadWrite( use_device );
    Execution::forAll( domain.numberOfDofs(), [=] MFEM_HOST_DEVICE( int dof ) {
      gap[dof] = area[dof] > 0.0 ? weighted_gap[dof] / area[dof] : 0.0;
    } );
    if ( state.dual != nullptr ) {
      if ( state.dual->Size() != domain.numberOfDofs() ) {
        throw std::invalid_argument( "SingleMortar multiplier has the wrong size." );
      }
      kinematics.pressure = *state.dual;
    }
  }

  void evaluateKinematics( const ContactDomain& domain, const InteractionStorage& interactions, const State& state,
                           Kinematics& kinematics ) const
  {
    evaluateKinematics<SequentialExecution>( domain, interactions, state, kinematics );
  }

  template <ExecutionPolicy Execution>
  Real addResidual( const ContactDomain& domain, const InteractionStorage& interactions, const State& state,
                    const Kinematics&, mfem::Vector& residual ) const
  {
    validateDual( domain, state );
    validateCoordinateVector( domain, residual );
    const bool use_device = Execution::uses_device;
    const auto view = domain.view( use_device );
    const auto* pairs = interactions.pairs.Read( use_device );
    const auto* dual = state.dual->Read( use_device );
    auto* output = residual.ReadWrite( use_device );
    const Parameters parameters = parameters_;
    forEachInteraction<Execution>( interactions.batches, interactions.coloring, [=] MFEM_HOST_DEVICE( int pair_index ) {
      using AD = detail::Gradient<Real, max_local_coordinates>;
      std::array<AD, max_local_coordinates> coordinates{};
      loadCoordinates( view, pairs[pair_index], coordinates, nullptr );
      const AD functional = multiplierFunctional( view, pairs[pair_index], coordinates, dual, parameters );
      scatterGradient( view, pairs[pair_index], functional, output );
    } );
    return 0.0;
  }

  Real addResidual( const ContactDomain& domain, const InteractionStorage& interactions, const State& state,
                    const Kinematics& kinematics, mfem::Vector& residual ) const
  {
    return addResidual<SequentialExecution>( domain, interactions, state, kinematics, residual );
  }

  template <ExecutionPolicy Execution>
  void applyCoordinateDerivative( const ContactDomain& domain, const InteractionStorage& interactions,
                                  const State& state, const Kinematics&, const mfem::Vector& direction,
                                  mfem::Vector& derivative ) const
  {
    validateDual( domain, state );
    validateCoordinateVector( domain, direction );
    validateCoordinateVector( domain, derivative );
    derivative = 0.0;
    const bool use_device = Execution::uses_device;
    const auto view = domain.view( use_device );
    const auto* pairs = interactions.pairs.Read( use_device );
    const auto* dual = state.dual->Read( use_device );
    const auto* vector = direction.Read( use_device );
    auto* output = derivative.ReadWrite( use_device );
    const Parameters parameters = parameters_;
    forEachInteraction<Execution>( interactions.batches, interactions.coloring, [=] MFEM_HOST_DEVICE( int pair_index ) {
      using Directional = detail::Gradient<Real, 1>;
      using AD = detail::Gradient<Directional, max_local_coordinates>;
      std::array<AD, max_local_coordinates> coordinates{};
      loadNestedCoordinates( view, pairs[pair_index], coordinates, vector );
      const AD functional = multiplierFunctional( view, pairs[pair_index], coordinates, dual, parameters );
      scatterHessianVector( view, pairs[pair_index], functional, output );
    } );
  }

  void applyCoordinateDerivative( const ContactDomain& domain, const InteractionStorage& interactions,
                                  const State& state, const Kinematics& kinematics, const mfem::Vector& direction,
                                  mfem::Vector& derivative ) const
  {
    applyCoordinateDerivative<SequentialExecution>( domain, interactions, state, kinematics, direction, derivative );
  }

  template <ExecutionPolicy Execution>
  void applyDualDerivative( const ContactDomain& domain, const InteractionStorage& interactions,
                            const mfem::Vector& dual_direction, mfem::Vector& force_direction ) const
  {
    if ( dual_direction.Size() != domain.numberOfDofs() ) {
      throw std::invalid_argument( "SingleMortar multiplier direction has the wrong size." );
    }
    validateCoordinateVector( domain, force_direction );
    force_direction = 0.0;
    const bool use_device = Execution::uses_device;
    const auto view = domain.view( use_device );
    const auto* pairs = interactions.pairs.Read( use_device );
    const auto* dual = dual_direction.Read( use_device );
    auto* output = force_direction.ReadWrite( use_device );
    const Parameters parameters = parameters_;
    forEachInteraction<Execution>( interactions.batches, interactions.coloring, [=] MFEM_HOST_DEVICE( int pair_index ) {
      using AD = detail::Gradient<Real, max_local_coordinates>;
      std::array<AD, max_local_coordinates> coordinates{};
      loadCoordinates( view, pairs[pair_index], coordinates, nullptr );
      const AD functional = multiplierFunctional( view, pairs[pair_index], coordinates, dual, parameters );
      scatterGradient( view, pairs[pair_index], functional, output );
    } );
  }

  void applyDualDerivative( const ContactDomain& domain, const InteractionStorage& interactions,
                            const mfem::Vector& dual_direction, mfem::Vector& force_direction ) const
  {
    applyDualDerivative<SequentialExecution>( domain, interactions, dual_direction, force_direction );
  }

  template <ExecutionPolicy Execution>
  void applyGapDerivative( const ContactDomain& domain, const InteractionStorage& interactions,
                           const Kinematics& kinematics, const mfem::Vector& direction,
                           mfem::Vector& gap_direction ) const
  {
    validateCoordinateVector( domain, direction );
    if ( gap_direction.Size() != domain.numberOfDofs() ) {
      throw std::invalid_argument( "SingleMortar gap derivative has the wrong size." );
    }
    const bool use_device = Execution::uses_device;
    gap_direction = 0.0;
    kinematics.directional_area = 0.0;
    gap_direction.UseDevice( use_device );
    kinematics.directional_area.UseDevice( use_device );
    auto* weighted_dot = gap_direction.ReadWrite( use_device );
    auto* area_dot = kinematics.directional_area.ReadWrite( use_device );
    const auto view = domain.view( use_device );
    const auto* pairs = interactions.pairs.Read( use_device );
    const auto* vector = direction.Read( use_device );
    const Parameters parameters = parameters_;
    forEachInteraction<Execution>( interactions.batches, interactions.coloring, [=] MFEM_HOST_DEVICE( int pair_index ) {
      using Directional = detail::Gradient<Real, 1>;
      std::array<Directional, max_local_coordinates> coordinates{};
      loadCoordinates( view, pairs[pair_index], coordinates, vector );
      const auto geometry = pairValues( view, pairs[pair_index], &coordinates, parameters );
      const int mortar = pairs[pair_index].mortar_element;
      const int nodes = std::min( view.numberOfElementDofs( mortar ), max_nodes_per_element );
      for ( int node = 0; node < nodes; ++node ) {
        const int dof = view.elementDof( mortar, node );
        AtomicAdd( weighted_dot[dof], geometry.weighted_gap[node].derivative[0] );
        AtomicAdd( area_dot[dof], geometry.area[node].derivative[0] );
      }
    } );
  }

  void applyGapDerivative( const ContactDomain& domain, const InteractionStorage& interactions,
                           const Kinematics& kinematics, const mfem::Vector& direction,
                           mfem::Vector& gap_direction ) const
  {
    applyGapDerivative<SequentialExecution>( domain, interactions, kinematics, direction, gap_direction );
  }

 private:
  template <typename T>
  struct PairValues {
    std::array<T, 3> mortar_centroid{};
    std::array<T, 3> nonmortar_centroid{};
    std::array<T, 3> normal{};
    std::array<T, max_nodes_per_element> weighted_gap{};
    std::array<T, max_nodes_per_element> area{};
    T gap{};
    T measure{};
  };

  template <typename T>
  struct ProjectedPolygon {
    std::array<std::array<T, 2>, 8> points{};
    int size{};
  };

  template <typename T>
  struct FaceValues {
    std::array<T, 3> point{};
    std::array<T, max_nodes_per_element> shape{};
  };

  static bool isFace( ElementTopology topology )
  {
    return topology == ElementTopology::Triangle || topology == ElementTopology::Quadrilateral;
  }

  template <typename T>
  TRIBOL_FUTURE_HOST_DEVICE static T squareRoot( const T& value )
  {
    using detail::sqrt;
    using std::sqrt;
    return sqrt( value );
  }

  template <typename T>
  TRIBOL_FUTURE_HOST_DEVICE static T cross2d( const std::array<T, 2>& left, const std::array<T, 2>& right )
  {
    return left[0] * right[1] - left[1] * right[0];
  }

  template <typename T>
  TRIBOL_FUTURE_HOST_DEVICE static std::array<T, 2> subtract2d( const std::array<T, 2>& left,
                                                                const std::array<T, 2>& right )
  {
    return { left[0] - right[0], left[1] - right[1] };
  }

  template <typename T>
  TRIBOL_FUTURE_HOST_DEVICE static std::array<T, 3> pointAt(
      const ContactDomainView& domain, int element, int side, int node,
      const std::array<T, max_local_coordinates>* coordinate_values )
  {
    std::array<T, 3> point{};
    const int dof = domain.elementDof( element, node );
    for ( int component = 0; component < 3; ++component ) {
      const int local = ( side * max_nodes_per_element + node ) * 3 + component;
      point[component] =
          coordinate_values == nullptr ? T{ domain.coordinate( component, dof ) } : ( *coordinate_values )[local];
    }
    return point;
  }

  template <typename T>
  TRIBOL_FUTURE_HOST_DEVICE static std::array<T, 3> faceNormal(
      const ContactDomainView& domain, int element, int side,
      const std::array<T, max_local_coordinates>* coordinate_values )
  {
    const int count = std::min( domain.numberOfElementDofs( element ), max_nodes_per_element );
    const auto first = subtract( pointAt( domain, element, side, 1, coordinate_values ),
                                 pointAt( domain, element, side, 0, coordinate_values ) );
    const auto second = subtract( pointAt( domain, element, side, 2, coordinate_values ),
                                  pointAt( domain, element, side, 0, coordinate_values ) );
    auto normal = cross( first, second );
    if ( count == 4 ) {
      const auto third = subtract( pointAt( domain, element, side, 3, coordinate_values ),
                                   pointAt( domain, element, side, 0, coordinate_values ) );
      const auto second_normal = cross( second, third );
      for ( int component = 0; component < 3; ++component ) {
        normal[component] = normal[component] + second_normal[component];
      }
    }
    const T magnitude = norm( normal );
    if ( std::abs( detail::scalarValue( magnitude ) ) <= 1.0e-28 ) {
      return {};
    }
    for ( int component = 0; component < 3; ++component ) {
      normal[component] = normal[component] / magnitude;
    }
    return normal;
  }

  template <typename T>
  TRIBOL_FUTURE_HOST_DEVICE static std::array<T, 2> projectToPlane( const std::array<T, 3>& point,
                                                                    const std::array<T, 3>& origin,
                                                                    const std::array<T, 3>& first_axis,
                                                                    const std::array<T, 3>& second_axis )
  {
    const auto difference = subtract( point, origin );
    return { difference[0] * first_axis[0] + difference[1] * first_axis[1] + difference[2] * first_axis[2],
             difference[0] * second_axis[0] + difference[1] * second_axis[1] + difference[2] * second_axis[2] };
  }

  template <typename T>
  TRIBOL_FUTURE_HOST_DEVICE static bool insideHalfPlane( const std::array<T, 2>& point, const std::array<T, 2>& first,
                                                         const std::array<T, 2>& second, Real orientation )
  {
    return orientation * detail::scalarValue( cross2d( subtract2d( second, first ), subtract2d( point, first ) ) ) >=
           -1.0e-14;
  }

  template <typename T>
  TRIBOL_FUTURE_HOST_DEVICE static std::array<T, 2> lineIntersection( const std::array<T, 2>& first,
                                                                      const std::array<T, 2>& second,
                                                                      const std::array<T, 2>& clip_first,
                                                                      const std::array<T, 2>& clip_second )
  {
    const auto direction = subtract2d( second, first );
    const auto clip_direction = subtract2d( clip_second, clip_first );
    const T denominator = cross2d( direction, clip_direction );
    if ( std::abs( detail::scalarValue( denominator ) ) < 1.0e-28 ) {
      return first;
    }
    const T parameter = cross2d( subtract2d( clip_first, first ), clip_direction ) / denominator;
    return { first[0] + parameter * direction[0], first[1] + parameter * direction[1] };
  }

  template <typename T>
  TRIBOL_FUTURE_HOST_DEVICE static ProjectedPolygon<T> clipPolygon(
      ProjectedPolygon<T> subject, const std::array<std::array<T, 2>, max_nodes_per_element>& clip, int clip_size )
  {
    T signed_area{};
    for ( int node = 0; node < clip_size; ++node ) {
      signed_area = signed_area + cross2d( clip[node], clip[( node + 1 ) % clip_size] );
    }
    const Real orientation = detail::scalarValue( signed_area ) >= 0.0 ? 1.0 : -1.0;
    for ( int edge = 0; edge < clip_size && subject.size > 0; ++edge ) {
      ProjectedPolygon<T> output;
      const auto clip_first = clip[edge];
      const auto clip_second = clip[( edge + 1 ) % clip_size];
      auto previous = subject.points[subject.size - 1];
      bool previous_inside = insideHalfPlane( previous, clip_first, clip_second, orientation );
      for ( int vertex = 0; vertex < subject.size; ++vertex ) {
        const auto current = subject.points[vertex];
        const bool current_inside = insideHalfPlane( current, clip_first, clip_second, orientation );
        if ( current_inside != previous_inside && output.size < static_cast<int>( output.points.size() ) ) {
          output.points[output.size++] = lineIntersection( previous, current, clip_first, clip_second );
        }
        if ( current_inside && output.size < static_cast<int>( output.points.size() ) ) {
          output.points[output.size++] = current;
        }
        previous = current;
        previous_inside = current_inside;
      }
      subject = output;
    }
    return subject;
  }

  template <typename T>
  TRIBOL_FUTURE_HOST_DEVICE static FaceValues<T> interpolateFace(
      const ContactDomainView& domain, int element, int side,
      const std::array<T, max_local_coordinates>* coordinate_values,
      const std::array<std::array<T, 2>, max_nodes_per_element>& projected_nodes,
      const std::array<T, 2>& projected_point )
  {
    FaceValues<T> result;
    const int count = std::min( domain.numberOfElementDofs( element ), max_nodes_per_element );
    if ( count == 3 ) {
      const auto first = subtract2d( projected_nodes[1], projected_nodes[0] );
      const auto second = subtract2d( projected_nodes[2], projected_nodes[0] );
      const auto relative = subtract2d( projected_point, projected_nodes[0] );
      const T denominator = cross2d( first, second );
      result.shape[1] = cross2d( relative, second ) / denominator;
      result.shape[2] = cross2d( first, relative ) / denominator;
      result.shape[0] = T{ 1.0 } - result.shape[1] - result.shape[2];
    } else {
      T reference_x{ 0.5 };
      T reference_y{ 0.5 };
      for ( int iteration = 0; iteration < 6; ++iteration ) {
        std::array<T, max_nodes_per_element> shape{};
        std::array<T, max_nodes_per_element> derivative_x{};
        std::array<T, max_nodes_per_element> derivative_y{};
        std::array<T, 2> mapped{};
        std::array<T, 2> tangent_x{};
        std::array<T, 2> tangent_y{};
        for ( int node = 0; node < count; ++node ) {
          const bool high_x = domain.referenceCoordinate( element, node, 0 ) > 0.5;
          const bool high_y = domain.referenceCoordinate( element, node, 1 ) > 0.5;
          const T factor_x = high_x ? reference_x : T{ 1.0 } - reference_x;
          const T factor_y = high_y ? reference_y : T{ 1.0 } - reference_y;
          shape[node] = factor_x * factor_y;
          derivative_x[node] = ( high_x ? T{ 1.0 } : T{ -1.0 } ) * factor_y;
          derivative_y[node] = factor_x * ( high_y ? T{ 1.0 } : T{ -1.0 } );
          for ( int component = 0; component < 2; ++component ) {
            mapped[component] = mapped[component] + shape[node] * projected_nodes[node][component];
            tangent_x[component] = tangent_x[component] + derivative_x[node] * projected_nodes[node][component];
            tangent_y[component] = tangent_y[component] + derivative_y[node] * projected_nodes[node][component];
          }
        }
        const auto residual = subtract2d( mapped, projected_point );
        const T determinant = cross2d( tangent_x, tangent_y );
        if ( std::abs( detail::scalarValue( determinant ) ) > 1.0e-28 ) {
          reference_x = reference_x - ( residual[0] * tangent_y[1] - residual[1] * tangent_y[0] ) / determinant;
          reference_y = reference_y - ( tangent_x[0] * residual[1] - tangent_x[1] * residual[0] ) / determinant;
        }
      }
      for ( int node = 0; node < count; ++node ) {
        const bool high_x = domain.referenceCoordinate( element, node, 0 ) > 0.5;
        const bool high_y = domain.referenceCoordinate( element, node, 1 ) > 0.5;
        result.shape[node] =
            ( high_x ? reference_x : T{ 1.0 } - reference_x ) * ( high_y ? reference_y : T{ 1.0 } - reference_y );
      }
    }
    for ( int node = 0; node < count; ++node ) {
      const auto point = pointAt( domain, element, side, node, coordinate_values );
      for ( int component = 0; component < 3; ++component ) {
        result.point[component] = result.point[component] + result.shape[node] * point[component];
      }
    }
    return result;
  }

  template <typename T>
  TRIBOL_FUTURE_HOST_DEVICE static T mortarBasis( const ContactDomainView& domain, int element, int node,
                                                  const std::array<T, max_nodes_per_element>& primal,
                                                  const BasisPolicy& basis )
  {
    const int count = domain.numberOfElementDofs( element );
    const bool high_x = domain.referenceCoordinate( element, node, 0 ) > 0.5;
    const bool high_y = domain.referenceCoordinate( element, node, 1 ) > 0.5;
    T reference_x{};
    T reference_y{};
    for ( int i = 0; i < count; ++i ) {
      reference_x = reference_x + T{ domain.referenceCoordinate( element, i, 0 ) } * primal[i];
      reference_y = reference_y + T{ domain.referenceCoordinate( element, i, 1 ) } * primal[i];
    }
    return basis.value( count, high_x, high_y, reference_x, reference_y, primal[node] );
  }

  template <typename T>
  TRIBOL_FUTURE_HOST_DEVICE static PairValues<T> pairValues(
      const ContactDomainView& domain, ElementPair pair, const std::array<T, max_local_coordinates>* coordinate_values,
      const Parameters& parameters )
  {
    PairValues<T> result;
    const int elements[2] = { pair.mortar_element, pair.nonmortar_element };
    const int counts[2] = { std::min( domain.numberOfElementDofs( elements[0] ), max_nodes_per_element ),
                            std::min( domain.numberOfElementDofs( elements[1] ), max_nodes_per_element ) };
    std::array<std::array<std::array<T, 3>, max_nodes_per_element>, 2> points{};
    for ( int side = 0; side < 2; ++side ) {
      for ( int node = 0; node < counts[side]; ++node ) {
        const int dof = domain.elementDof( elements[side], node );
        for ( int component = 0; component < 3; ++component ) {
          const int local = ( side * max_nodes_per_element + node ) * 3 + component;
          points[side][node][component] =
              coordinate_values == nullptr ? T{ domain.coordinate( component, dof ) } : ( *coordinate_values )[local];
        }
      }
    }

    result.normal = faceNormal( domain, elements[1], 1, coordinate_values );
    if ( std::abs( detail::scalarValue( norm( result.normal ) ) ) <= 1.0e-28 ) {
      return result;
    }
    auto first_axis = subtract( points[1][1], points[1][0] );
    const T first_axis_length = norm( first_axis );
    if ( std::abs( detail::scalarValue( first_axis_length ) ) <= 1.0e-28 ) {
      return result;
    }
    for ( int component = 0; component < 3; ++component ) {
      first_axis[component] = first_axis[component] / first_axis_length;
    }
    const auto second_axis = cross( result.normal, first_axis );
    const auto origin = points[1][0];
    std::array<std::array<std::array<T, 2>, max_nodes_per_element>, 2> projected{};
    for ( int side = 0; side < 2; ++side ) {
      for ( int node = 0; node < counts[side]; ++node ) {
        projected[side][node] = projectToPlane( points[side][node], origin, first_axis, second_axis );
      }
    }
    ProjectedPolygon<T> overlap;
    overlap.size = counts[0];
    for ( int node = 0; node < counts[0]; ++node ) {
      overlap.points[node] = projected[0][node];
    }
    overlap = clipPolygon( overlap, projected[1], counts[1] );
    if ( overlap.size < 3 ) {
      return result;
    }
    std::array<T, 2> polygon_centroid{};
    for ( int vertex = 0; vertex < overlap.size; ++vertex ) {
      polygon_centroid[0] = polygon_centroid[0] + overlap.points[vertex][0] / T{ Real( overlap.size ) };
      polygon_centroid[1] = polygon_centroid[1] + overlap.points[vertex][1] / T{ Real( overlap.size ) };
    }
    for ( int vertex = 0; vertex < overlap.size; ++vertex ) {
      const auto first = subtract2d( overlap.points[vertex], polygon_centroid );
      const auto second = subtract2d( overlap.points[( vertex + 1 ) % overlap.size], polygon_centroid );
      T triangle_area = T{ 0.5 } * cross2d( first, second );
      if ( detail::scalarValue( triangle_area ) < 0.0 ) {
        triangle_area = -triangle_area;
      }
      const int number_of_points = parameters.quadrature_points == 1 ? 1 : 3;
      for ( int point = 0; point < number_of_points; ++point ) {
        Real first_coordinate = 1.0 / 3.0;
        Real second_coordinate = 1.0 / 3.0;
        if ( number_of_points == 3 ) {
          first_coordinate = point == 1 ? 2.0 / 3.0 : 1.0 / 6.0;
          second_coordinate = point == 2 ? 2.0 / 3.0 : 1.0 / 6.0;
        }
        const Real center_coordinate = 1.0 - first_coordinate - second_coordinate;
        const std::array<T, 2> quadrature_point{
            T{ center_coordinate } * polygon_centroid[0] + T{ first_coordinate } * overlap.points[vertex][0] +
                T{ second_coordinate } * overlap.points[( vertex + 1 ) % overlap.size][0],
            T{ center_coordinate } * polygon_centroid[1] + T{ first_coordinate } * overlap.points[vertex][1] +
                T{ second_coordinate } * overlap.points[( vertex + 1 ) % overlap.size][1] };
        const T weight = triangle_area / T{ Real( number_of_points ) };
        const auto mortar =
            interpolateFace( domain, elements[0], 0, coordinate_values, projected[0], quadrature_point );
        const auto nonmortar =
            interpolateFace( domain, elements[1], 1, coordinate_values, projected[1], quadrature_point );
        T gap{};
        for ( int component = 0; component < 3; ++component ) {
          gap = gap + ( mortar.point[component] - nonmortar.point[component] ) * result.normal[component];
          result.mortar_centroid[component] = result.mortar_centroid[component] + weight * mortar.point[component];
          result.nonmortar_centroid[component] =
              result.nonmortar_centroid[component] + weight * nonmortar.point[component];
        }
        result.gap = result.gap + weight * gap;
        result.measure = result.measure + weight;
        for ( int node = 0; node < counts[0]; ++node ) {
          const T basis = mortarBasis( domain, elements[0], node, mortar.shape, parameters.basis );
          result.weighted_gap[node] = result.weighted_gap[node] + weight * basis * gap;
          result.area[node] = result.area[node] + weight * basis;
        }
      }
    }
    if ( std::abs( detail::scalarValue( result.measure ) ) > 1.0e-28 ) {
      result.gap = result.gap / result.measure;
      for ( int component = 0; component < 3; ++component ) {
        result.mortar_centroid[component] = result.mortar_centroid[component] / result.measure;
        result.nonmortar_centroid[component] = result.nonmortar_centroid[component] / result.measure;
      }
    }
    return result;
  }

  template <typename T>
  TRIBOL_FUTURE_HOST_DEVICE static T multiplierFunctional( const ContactDomainView& domain, ElementPair pair,
                                                           const std::array<T, max_local_coordinates>& coordinates,
                                                           const Real* dual, const Parameters& parameters )
  {
    const auto geometry = pairValues( domain, pair, &coordinates, parameters );
    const int mortar = pair.mortar_element;
    const int nodes = std::min( domain.numberOfElementDofs( mortar ), max_nodes_per_element );
    T result{};
    for ( int node = 0; node < nodes; ++node ) {
      result = result + T{ dual[domain.elementDof( mortar, node )] } * geometry.weighted_gap[node];
    }
    return result;
  }

  template <typename T>
  TRIBOL_FUTURE_HOST_DEVICE static std::array<T, 3> subtract( const std::array<T, 3>& left,
                                                              const std::array<T, 3>& right )
  {
    return { left[0] - right[0], left[1] - right[1], left[2] - right[2] };
  }

  template <typename T>
  TRIBOL_FUTURE_HOST_DEVICE static std::array<T, 3> cross( const std::array<T, 3>& left, const std::array<T, 3>& right )
  {
    return { left[1] * right[2] - left[2] * right[1], left[2] * right[0] - left[0] * right[2],
             left[0] * right[1] - left[1] * right[0] };
  }

  template <typename T>
  TRIBOL_FUTURE_HOST_DEVICE static T norm( const std::array<T, 3>& vector )
  {
    return squareRoot( vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2] );
  }

  template <typename T>
  TRIBOL_FUTURE_HOST_DEVICE static void loadCoordinates( const ContactDomainView& domain, ElementPair pair,
                                                         std::array<T, max_local_coordinates>& coordinates,
                                                         const Real* direction )
  {
    const int elements[2] = { pair.mortar_element, pair.nonmortar_element };
    for ( int side = 0; side < 2; ++side ) {
      const int nodes = std::min( domain.numberOfElementDofs( elements[side] ), max_nodes_per_element );
      for ( int node = 0; node < nodes; ++node ) {
        const int dof = domain.elementDof( elements[side], node );
        for ( int component = 0; component < 3; ++component ) {
          const int local = ( side * max_nodes_per_element + node ) * 3 + component;
          const int global = component * domain.number_of_dofs + dof;
          if constexpr ( std::is_same_v<T, Real> ) {
            coordinates[local] = domain.coordinates[global];
          } else if ( direction == nullptr ) {
            coordinates[local] = T::variable( domain.coordinates[global], local );
          } else {
            coordinates[local] = T{ domain.coordinates[global] };
            coordinates[local].derivative[0] = direction[global];
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
      const int nodes = std::min( domain.numberOfElementDofs( elements[side] ), max_nodes_per_element );
      for ( int node = 0; node < nodes; ++node ) {
        const int dof = domain.elementDof( elements[side], node );
        for ( int component = 0; component < 3; ++component ) {
          const int local = ( side * max_nodes_per_element + node ) * 3 + component;
          const int global = component * domain.number_of_dofs + dof;
          Directional value( domain.coordinates[global] );
          value.derivative[0] = direction[global];
          coordinates[local] = T::variable( value, local );
        }
      }
    }
  }

  template <typename AD>
  TRIBOL_FUTURE_HOST_DEVICE static void scatterGradient( const ContactDomainView& domain, ElementPair pair,
                                                         const AD& functional, Real* output )
  {
    const int elements[2] = { pair.mortar_element, pair.nonmortar_element };
    for ( int side = 0; side < 2; ++side ) {
      const int nodes = std::min( domain.numberOfElementDofs( elements[side] ), max_nodes_per_element );
      for ( int node = 0; node < nodes; ++node ) {
        const int dof = domain.elementDof( elements[side], node );
        for ( int component = 0; component < 3; ++component ) {
          const int local = ( side * max_nodes_per_element + node ) * 3 + component;
          AtomicAdd( output[component * domain.number_of_dofs + dof],
                     detail::scalarValue( functional.derivative[local] ) );
        }
      }
    }
  }

  template <typename AD>
  TRIBOL_FUTURE_HOST_DEVICE static void scatterHessianVector( const ContactDomainView& domain, ElementPair pair,
                                                              const AD& functional, Real* output )
  {
    const int elements[2] = { pair.mortar_element, pair.nonmortar_element };
    for ( int side = 0; side < 2; ++side ) {
      const int nodes = std::min( domain.numberOfElementDofs( elements[side] ), max_nodes_per_element );
      for ( int node = 0; node < nodes; ++node ) {
        const int dof = domain.elementDof( elements[side], node );
        for ( int component = 0; component < 3; ++component ) {
          const int local = ( side * max_nodes_per_element + node ) * 3 + component;
          AtomicAdd( output[component * domain.number_of_dofs + dof], functional.derivative[local].derivative[0] );
        }
      }
    }
  }

  static void validateDual( const ContactDomain& domain, const State& state )
  {
    if ( state.dual == nullptr || state.dual->Size() != domain.numberOfDofs() ) {
      throw std::invalid_argument( "SingleMortar requires a contact-space multiplier vector." );
    }
  }

  static void validateCoordinateVector( const ContactDomain& domain, const mfem::Vector& vector )
  {
    if ( vector.Size() != domain.coordinateSize() ) {
      throw std::invalid_argument( "SingleMortar coordinate vector has the wrong size." );
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
    kinematics.directional_area.UseDevice( use_device );
    const_cast<mfem::Array<ElementPair>&>( interactions.pairs ).GetMemory().UseDevice( use_device );
    const_cast<InteractionColoring&>( interactions.coloring ).offsets.GetMemory().UseDevice( use_device );
    const_cast<InteractionColoring&>( interactions.coloring ).pair_indices.GetMemory().UseDevice( use_device );
  }

  Parameters parameters_;
};

}  // namespace tribol::future

#endif
