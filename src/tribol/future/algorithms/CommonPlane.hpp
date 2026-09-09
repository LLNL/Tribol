#ifndef TRIBOL_FUTURE_ALGORITHMS_COMMONPLANE_HPP_
#define TRIBOL_FUTURE_ALGORITHMS_COMMONPLANE_HPP_

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
#include <type_traits>
#include <vector>

namespace tribol::future {

template <CommonPlaneNormalLaw NormalLaw = ConstantKinematicPenalty, CommonPlaneRateLaw RateLaw = NoRatePenalty,
          CommonPlaneTangentialLaw TangentialLaw = FrictionlessTangentialResponse,
          LinearizationPolicy LinearizationPolicy = AnalyticLinearization>
class CommonPlane {
 public:
  static constexpr int max_nodes_per_element = 4;
  static constexpr int max_local_coordinates = 2 * max_nodes_per_element * 3;
  static constexpr int maximum_element_order = 1;

  struct Parameters {
    NormalLaw normal{};
    RateLaw rate{};
    TangentialLaw tangential{};
    Real minimum_measure{ 1.0e-14 };
  };

  struct State {
    const mfem::Vector* velocity{};
    const mfem::Vector* reference_coordinates{};
    const mfem::Vector* element_thickness{};
    const mfem::Vector* material_modulus{};
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

  static constexpr bool accepts_external_pressure = false;
  static constexpr bool has_nodal_kinematics = false;
  static constexpr MethodCapabilities capabilities{
      .needs_dual_space = false,
      .needs_velocity = RateLaw::needs_velocity || TangentialLaw::needs_velocity,
      .needs_reference_field = NormalLaw::needs_reference_field,
      .needs_material_fields = NormalLaw::needs_material_fields,
      .has_history = false,
      .has_energy = RateLaw::is_conservative && TangentialLaw::is_conservative,
      .supports_self_contact = true,
      .supports_assembled_jacobian = true,
  };

  explicit CommonPlane( Parameters parameters = {} ) : parameters_( parameters )
  {
    if ( parameters_.minimum_measure <= 0.0 ) {
      throw std::invalid_argument( "CommonPlane minimum measure must be positive." );
    }
  }

  void buildInteractions( const ContactDomain& domain, const CandidatePairs& candidates,
                          InteractionStorage& interactions ) const
  {
    const auto view = domain.view( false );
    std::vector<ElementPair> accepted;
    accepted.reserve( candidates.pairs.Size() );
    for ( int pair_index = 0; pair_index < candidates.pairs.Size(); ++pair_index ) {
      const auto pair = candidates.pairs[pair_index];
      if ( view.orders[pair.mortar_element] != 1 || view.orders[pair.nonmortar_element] != 1 ) {
        throw std::invalid_argument( "CommonPlane requires linear surface elements; use LowOrderRefinedSurface." );
      }
      const bool supported =
          view.dimension == 2
              ? view.topology( pair.mortar_element ) == ElementTopology::Segment &&
                    view.topology( pair.nonmortar_element ) == ElementTopology::Segment
              : isFace( view.topology( pair.mortar_element ) ) && isFace( view.topology( pair.nonmortar_element ) );
      if ( supported && pairGeometry<Real>( view, pair ).measure > parameters_.minimum_measure ) {
        accepted.push_back( pair );
      }
    }
    buildInteractionBatches( domain, std::move( accepted ), interactions.pairs, interactions.batches );
    buildInteractionColoring( domain, interactions.pairs, interactions.coloring );
  }

  void prepare( const ContactDomain& domain, const InteractionStorage& interactions, Kinematics& kinematics ) const
  {
    const int pairs = interactions.pairs.Size();
    kinematics.gap.SetSize( pairs );
    kinematics.weighted_gap.SetSize( 0 );
    kinematics.tributary_area.SetSize( 0 );
    kinematics.pressure.SetSize( pairs );
    kinematics.pressure_tangent.SetSize( pairs );
    kinematics.potential_density.SetSize( pairs );
    kinematics.pair_offsets.SetSize( pairs + 1 );
    kinematics.quadrature_points.SetSize( pairs * 2 * domain.dimension() );
    kinematics.quadrature_weights.SetSize( pairs );
    kinematics.quadrature_gaps.SetSize( pairs );
    kinematics.scalar_workspace.SetSize( interactions.pairs.Size() + 1 );
    for ( int pair = 0; pair <= pairs; ++pair ) {
      kinematics.pair_offsets[pair] = pair;
    }
  }

  template <ExecutionPolicy Execution>
  void evaluateKinematics( const ContactDomain& domain, const InteractionStorage& interactions, const State&,
                           Kinematics& kinematics ) const
  {
    const bool use_device = Execution::uses_device;
    setMemoryMode( kinematics, interactions, use_device );
    const auto view = domain.view( use_device );
    const auto* pairs = interactions.pairs.Read( use_device );
    auto* gap = kinematics.gap.Write( use_device );
    auto* weights = kinematics.quadrature_weights.Write( use_device );
    auto* gaps = kinematics.quadrature_gaps.Write( use_device );
    auto* points = kinematics.quadrature_points.Write( use_device );
    const int dimension = domain.dimension();
    forEachInteraction<Execution>( interactions.batches, interactions.coloring, [=] MFEM_HOST_DEVICE( int pair_index ) {
      const auto geometry = pairGeometry<Real>( view, pairs[pair_index] );
      gap[pair_index] = geometry.gap;
      gaps[pair_index] = geometry.gap;
      weights[pair_index] = geometry.measure;
      for ( int component = 0; component < dimension; ++component ) {
        points[( 2 * pair_index ) * dimension + component] = geometry.first_centroid[component];
        points[( 2 * pair_index + 1 ) * dimension + component] = geometry.second_centroid[component];
      }
    } );
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
    validateState( domain, state );
    validateCoordinateVector( domain, residual );
    const bool use_device = Execution::uses_device;
    const auto view = domain.view( use_device );
    const auto* pairs = interactions.pairs.Read( use_device );
    const Real* velocity = state.velocity != nullptr ? state.velocity->Read( use_device ) : nullptr;
    const Real* reference_coordinates =
        state.reference_coordinates != nullptr ? state.reference_coordinates->Read( use_device ) : nullptr;
    const Real* thickness = state.element_thickness != nullptr ? state.element_thickness->Read( use_device ) : nullptr;
    const Real* modulus = state.material_modulus != nullptr ? state.material_modulus->Read( use_device ) : nullptr;
    Real* output = residual.ReadWrite( use_device );
    auto& energy_workspace = kinematics.scalar_workspace;
    energy_workspace.UseDevice( use_device );
    energy_workspace = 0.0;
    Real* energy = energy_workspace.ReadWrite( use_device );
    const Parameters parameters = parameters_;
    forEachInteraction<Execution>( interactions.batches, interactions.coloring, [=] MFEM_HOST_DEVICE( int pair_index ) {
      using AD = detail::Gradient<Real, max_local_coordinates>;
      std::array<AD, max_local_coordinates> coordinates{};
      const int local_size = loadCoordinates( view, pairs[pair_index], coordinates, nullptr );
      const AD pair_energy = contactPotential( view, pairs[pair_index], coordinates, parameters, reference_coordinates,
                                               thickness, modulus );
      scatterGradient( view, pairs[pair_index], pair_energy, local_size, output );
      addVelocityResponse( view, pairs[pair_index], parameters, velocity, nullptr, thickness, modulus, output );
      if constexpr ( capabilities.has_energy ) {
        storeInteractionScalar<Execution>( energy, pair_index, detail::scalarValue( pair_energy.value ) );
      }
    } );
    if constexpr ( capabilities.has_energy ) {
      finalizeInteractionScalar<Execution>( kinematics.scalar_workspace, interactions.pairs.Size() );
      Execution::synchronize();
      return energy_workspace.HostRead()[0];
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
                                  const State& state, const Kinematics&, const mfem::Vector& direction,
                                  mfem::Vector& derivative ) const
  {
    validateState( domain, state );
    validateCoordinateVector( domain, direction );
    validateCoordinateVector( domain, derivative );
    derivative = 0.0;
    const bool use_device = Execution::uses_device;
    const auto view = domain.view( use_device );
    const auto* pairs = interactions.pairs.Read( use_device );
    const auto* vector = direction.Read( use_device );
    const Real* velocity = state.velocity != nullptr ? state.velocity->Read( use_device ) : nullptr;
    const Real* reference_coordinates =
        state.reference_coordinates != nullptr ? state.reference_coordinates->Read( use_device ) : nullptr;
    const Real* thickness = state.element_thickness != nullptr ? state.element_thickness->Read( use_device ) : nullptr;
    const Real* modulus = state.material_modulus != nullptr ? state.material_modulus->Read( use_device ) : nullptr;
    auto* output = derivative.ReadWrite( use_device );
    const Parameters parameters = parameters_;
    forEachInteraction<Execution>( interactions.batches, interactions.coloring, [=] MFEM_HOST_DEVICE( int pair_index ) {
      using Directional = detail::Gradient<Real, 1>;
      using AD = detail::Gradient<Directional, max_local_coordinates>;
      std::array<AD, max_local_coordinates> coordinates{};
      const int local_size = loadNestedCoordinates( view, pairs[pair_index], coordinates, vector );
      const AD pair_energy = contactPotential( view, pairs[pair_index], coordinates, parameters, reference_coordinates,
                                               thickness, modulus );
      scatterHessianVector( view, pairs[pair_index], pair_energy, local_size, output );
      if constexpr ( capabilities.needs_velocity ) {
        std::array<Directional, max_local_coordinates> directional_coordinates{};
        loadCoordinates( view, pairs[pair_index], directional_coordinates, vector );
        addVelocityCoordinateDerivative( view, pairs[pair_index], parameters, velocity, thickness, modulus,
                                         directional_coordinates, output );
      }
    } );
  }

  void applyCoordinateDerivative( const ContactDomain& domain, const InteractionStorage& interactions,
                                  const State& state, const Kinematics& kinematics, const mfem::Vector& direction,
                                  mfem::Vector& derivative ) const
  {
    applyCoordinateDerivative<SequentialExecution>( domain, interactions, state, kinematics, direction, derivative );
  }

  template <ExecutionPolicy Execution>
  void applyVelocityDerivative( const ContactDomain& domain, const InteractionStorage& interactions, const State& state,
                                const mfem::Vector& direction, mfem::Vector& derivative ) const
  {
    validateState( domain, state );
    validateCoordinateVector( domain, direction );
    validateCoordinateVector( domain, derivative );
    derivative = 0.0;
    const bool use_device = Execution::uses_device;
    const auto view = domain.view( use_device );
    const auto* pairs = interactions.pairs.Read( use_device );
    const Real* velocity = state.velocity != nullptr ? state.velocity->Read( use_device ) : nullptr;
    const Real* velocity_direction = direction.Read( use_device );
    const Real* thickness = state.element_thickness != nullptr ? state.element_thickness->Read( use_device ) : nullptr;
    const Real* modulus = state.material_modulus != nullptr ? state.material_modulus->Read( use_device ) : nullptr;
    Real* output = derivative.ReadWrite( use_device );
    const Parameters parameters = parameters_;
    forEachInteraction<Execution>( interactions.batches, interactions.coloring, [=] MFEM_HOST_DEVICE( int pair_index ) {
      addVelocityResponse( view, pairs[pair_index], parameters, velocity, velocity_direction, thickness, modulus,
                           output );
    } );
  }

  void applyVelocityDerivative( const ContactDomain& domain, const InteractionStorage& interactions, const State& state,
                                const mfem::Vector& direction, mfem::Vector& derivative ) const
  {
    applyVelocityDerivative<SequentialExecution>( domain, interactions, state, direction, derivative );
  }

 private:
  template <typename T>
  struct PairGeometry {
    std::array<T, 3> first_centroid{};
    std::array<T, 3> second_centroid{};
    std::array<T, 3> normal{};
    T gap{};
    T measure{};
  };

  template <typename T>
  struct ProjectedPolygon {
    std::array<std::array<T, 2>, 8> points{};
    int size{};
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
  TRIBOL_FUTURE_HOST_DEVICE static PairGeometry<T> pairGeometry(
      const ContactDomainView& domain, ElementPair pair, const std::array<T, max_local_coordinates>* x = nullptr )
  {
    PairGeometry<T> result;
    const int elements[2] = { pair.mortar_element, pair.nonmortar_element };
    const int counts[2] = { std::min( domain.numberOfElementDofs( elements[0] ), max_nodes_per_element ),
                            std::min( domain.numberOfElementDofs( elements[1] ), max_nodes_per_element ) };
    std::array<std::array<std::array<T, 3>, max_nodes_per_element>, 2> points{};
    for ( int side = 0; side < 2; ++side ) {
      for ( int node = 0; node < counts[side]; ++node ) {
        const int dof = domain.elementDof( elements[side], node );
        for ( int component = 0; component < domain.dimension; ++component ) {
          const int local = ( side * max_nodes_per_element + node ) * 3 + component;
          points[side][node][component] = x == nullptr ? T{ domain.coordinate( component, dof ) } : ( *x )[local];
        }
      }
    }

    std::array<std::array<T, 3>, 2> normals{};
    for ( int side = 0; side < 2; ++side ) {
      if ( domain.dimension == 2 ) {
        const T dx = points[side][1][0] - points[side][0][0];
        const T dy = points[side][1][1] - points[side][0][1];
        const T length = squareRoot( dx * dx + dy * dy );
        if ( detail::scalarValue( length ) <= 1.0e-28 ) {
          return result;
        }
        normals[side][0] = dy / length;
        normals[side][1] = -dx / length;
      } else {
        const auto first = subtract( points[side][1], points[side][0] );
        const auto second = subtract( points[side][2], points[side][0] );
        auto normal = cross( first, second );
        if ( counts[side] == 4 ) {
          const auto third = subtract( points[side][3], points[side][0] );
          const auto second_normal = cross( second, third );
          for ( int component = 0; component < 3; ++component ) {
            normal[component] = normal[component] + second_normal[component];
          }
        }
        const T normal_length = vectorNorm( normal );
        if ( detail::scalarValue( normal_length ) <= 1.0e-28 ) {
          return result;
        }
        for ( int component = 0; component < 3; ++component ) {
          normals[side][component] = normal[component] / normal_length;
        }
      }
    }
    std::array<T, 3> common_normal{};
    T normal_norm{};
    for ( int component = 0; component < domain.dimension; ++component ) {
      common_normal[component] = normals[1][component] - normals[0][component];
      normal_norm = normal_norm + common_normal[component] * common_normal[component];
    }
    normal_norm = squareRoot( normal_norm );
    if ( detail::scalarValue( normal_norm ) < 1.0e-14 ) {
      common_normal = normals[1];
    } else {
      for ( int component = 0; component < domain.dimension; ++component ) {
        common_normal[component] = common_normal[component] / normal_norm;
      }
    }
    result.normal = common_normal;

    if ( domain.dimension == 2 ) {
      const std::array<T, 2> tangent{ -common_normal[1], common_normal[0] };
      std::array<std::array<T, 2>, 2> projected{};
      for ( int side = 0; side < 2; ++side ) {
        for ( int node = 0; node < 2; ++node ) {
          projected[side][node] = points[side][node][0] * tangent[0] + points[side][node][1] * tangent[1];
        }
      }
      const T first_lower = detail::scalarValue( projected[0][0] ) < detail::scalarValue( projected[0][1] )
                                ? projected[0][0]
                                : projected[0][1];
      const T first_upper = detail::scalarValue( projected[0][0] ) < detail::scalarValue( projected[0][1] )
                                ? projected[0][1]
                                : projected[0][0];
      const T second_lower = detail::scalarValue( projected[1][0] ) < detail::scalarValue( projected[1][1] )
                                 ? projected[1][0]
                                 : projected[1][1];
      const T second_upper = detail::scalarValue( projected[1][0] ) < detail::scalarValue( projected[1][1] )
                                 ? projected[1][1]
                                 : projected[1][0];
      const T lower =
          detail::scalarValue( first_lower ) > detail::scalarValue( second_lower ) ? first_lower : second_lower;
      const T upper =
          detail::scalarValue( first_upper ) < detail::scalarValue( second_upper ) ? first_upper : second_upper;
      if ( detail::scalarValue( upper ) <= detail::scalarValue( lower ) ) {
        return result;
      }
      result.measure = upper - lower;
      const T midpoint = T{ 0.5 } * ( lower + upper );
      for ( int side = 0; side < 2; ++side ) {
        const T parameter = ( midpoint - projected[side][0] ) / ( projected[side][1] - projected[side][0] );
        auto& centroid = side == 0 ? result.first_centroid : result.second_centroid;
        for ( int component = 0; component < 2; ++component ) {
          centroid[component] =
              points[side][0][component] + parameter * ( points[side][1][component] - points[side][0][component] );
        }
      }
    } else {
      auto first_axis = subtract( points[0][1], points[0][0] );
      T normal_component{};
      for ( int component = 0; component < 3; ++component ) {
        normal_component = normal_component + first_axis[component] * common_normal[component];
      }
      for ( int component = 0; component < 3; ++component ) {
        first_axis[component] = first_axis[component] - normal_component * common_normal[component];
      }
      const T axis_length = vectorNorm( first_axis );
      if ( detail::scalarValue( axis_length ) <= 1.0e-28 ) {
        return result;
      }
      for ( int component = 0; component < 3; ++component ) {
        first_axis[component] = first_axis[component] / axis_length;
      }
      const auto second_axis = cross( common_normal, first_axis );
      const auto origin = points[1][0];
      std::array<std::array<std::array<T, 2>, max_nodes_per_element>, 2> projected{};
      for ( int side = 0; side < 2; ++side ) {
        for ( int node = 0; node < counts[side]; ++node ) {
          const auto relative = subtract( points[side][node], origin );
          for ( int component = 0; component < 3; ++component ) {
            projected[side][node][0] = projected[side][node][0] + relative[component] * first_axis[component];
            projected[side][node][1] = projected[side][node][1] + relative[component] * second_axis[component];
          }
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
      T twice_area{};
      std::array<T, 2> projected_centroid{};
      for ( int vertex = 0; vertex < overlap.size; ++vertex ) {
        const auto& first = overlap.points[vertex];
        const auto& second = overlap.points[( vertex + 1 ) % overlap.size];
        const T cross_value = cross2d( first, second );
        twice_area = twice_area + cross_value;
        projected_centroid[0] = projected_centroid[0] + ( first[0] + second[0] ) * cross_value;
        projected_centroid[1] = projected_centroid[1] + ( first[1] + second[1] ) * cross_value;
      }
      if ( std::abs( detail::scalarValue( twice_area ) ) <= 1.0e-28 ) {
        return result;
      }
      projected_centroid[0] = projected_centroid[0] / ( T{ 3.0 } * twice_area );
      projected_centroid[1] = projected_centroid[1] / ( T{ 3.0 } * twice_area );
      result.measure = detail::scalarValue( twice_area ) < 0.0 ? T{ -0.5 } * twice_area : T{ 0.5 } * twice_area;
      std::array<T, 3> projected_point{};
      for ( int component = 0; component < 3; ++component ) {
        projected_point[component] = origin[component] + projected_centroid[0] * first_axis[component] +
                                     projected_centroid[1] * second_axis[component];
      }
      for ( int side = 0; side < 2; ++side ) {
        auto& centroid = side == 0 ? result.first_centroid : result.second_centroid;
        T numerator{};
        T denominator{};
        for ( int component = 0; component < 3; ++component ) {
          numerator =
              numerator + normals[side][component] * ( points[side][0][component] - projected_point[component] );
          denominator = denominator + normals[side][component] * common_normal[component];
        }
        if ( std::abs( detail::scalarValue( denominator ) ) <= 1.0e-28 ) {
          return PairGeometry<T>{};
        }
        const T distance = numerator / denominator;
        for ( int component = 0; component < 3; ++component ) {
          centroid[component] = projected_point[component] + distance * common_normal[component];
        }
      }
    }

    for ( int component = 0; component < domain.dimension; ++component ) {
      result.gap = result.gap +
                   ( result.first_centroid[component] - result.second_centroid[component] ) * common_normal[component];
    }
    return result;
  }

  template <typename T>
  TRIBOL_FUTURE_HOST_DEVICE static T contactPotential( const ContactDomainView& domain, ElementPair pair,
                                                       const std::array<T, max_local_coordinates>& coordinates,
                                                       const Parameters& parameters, const Real* reference_coordinates,
                                                       const Real* thickness, const Real* modulus )
  {
    const auto geometry = pairGeometry( domain, pair, &coordinates );
    const Real stiffness = normalStiffness( parameters.normal, thickness, modulus, pair );
    if constexpr ( tiedNormal ) {
      const int elements[2] = { pair.mortar_element, pair.nonmortar_element };
      std::array<T, 3> displacement_jump{};
      for ( int side = 0; side < 2; ++side ) {
        const int count = std::min( domain.numberOfElementDofs( elements[side] ), max_nodes_per_element );
        const T sign = side == 0 ? T{ 1.0 } : T{ -1.0 };
        for ( int node = 0; node < count; ++node ) {
          const int dof = domain.elementDof( elements[side], node );
          for ( int component = 0; component < domain.dimension; ++component ) {
            const int local = ( side * max_nodes_per_element + node ) * 3 + component;
            const T displacement =
                coordinates[local] - T{ reference_coordinates[component * domain.number_of_dofs + dof] };
            displacement_jump[component] =
                displacement_jump[component] + sign * displacement / T{ static_cast<Real>( count ) };
          }
        }
      }
      T displacement_squared{};
      if constexpr ( NormalLaw::ties_tangential ) {
        for ( int component = 0; component < domain.dimension; ++component ) {
          displacement_squared = displacement_squared + displacement_jump[component] * displacement_jump[component];
        }
      } else {
        T normal_displacement{};
        for ( int component = 0; component < domain.dimension; ++component ) {
          normal_displacement = normal_displacement + displacement_jump[component] * geometry.normal[component];
        }
        displacement_squared = normal_displacement * normal_displacement;
      }
      return T{ 0.5 * stiffness } * geometry.measure * displacement_squared;
    } else {
      const T active_gap = detail::activeMinimum( geometry.gap );
      return T{ 0.5 * stiffness } * geometry.measure * active_gap * active_gap;
    }
  }

  TRIBOL_FUTURE_HOST_DEVICE static Real normalStiffness( const NormalLaw& law, const Real* thickness,
                                                         const Real* modulus, ElementPair pair )
  {
    return law.penaltyStiffness( pair.mortar_element, pair.nonmortar_element, thickness, modulus );
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
  TRIBOL_FUTURE_HOST_DEVICE static T vectorNorm( const std::array<T, 3>& vector )
  {
    return squareRoot( vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2] );
  }

  template <typename T>
  TRIBOL_FUTURE_HOST_DEVICE static int loadCoordinates( const ContactDomainView& domain, ElementPair pair,
                                                        std::array<T, max_local_coordinates>& coordinates,
                                                        const Real* direction )
  {
    const int elements[2] = { pair.mortar_element, pair.nonmortar_element };
    int local_size{};
    for ( int side = 0; side < 2; ++side ) {
      const int count = std::min( domain.numberOfElementDofs( elements[side] ), max_nodes_per_element );
      for ( int node = 0; node < count; ++node ) {
        const int dof = domain.elementDof( elements[side], node );
        for ( int component = 0; component < domain.dimension; ++component ) {
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
          local_size = std::max( local_size, local + 1 );
        }
      }
    }
    return local_size;
  }

  template <typename T>
  TRIBOL_FUTURE_HOST_DEVICE static int loadNestedCoordinates( const ContactDomainView& domain, ElementPair pair,
                                                              std::array<T, max_local_coordinates>& coordinates,
                                                              const Real* direction )
  {
    using Directional = decltype( T{}.value );
    const int elements[2] = { pair.mortar_element, pair.nonmortar_element };
    int local_size{};
    for ( int side = 0; side < 2; ++side ) {
      const int count = std::min( domain.numberOfElementDofs( elements[side] ), max_nodes_per_element );
      for ( int node = 0; node < count; ++node ) {
        const int dof = domain.elementDof( elements[side], node );
        for ( int component = 0; component < domain.dimension; ++component ) {
          const int local = ( side * max_nodes_per_element + node ) * 3 + component;
          const int global = component * domain.number_of_dofs + dof;
          Directional value( domain.coordinates[global] );
          value.derivative[0] = direction[global];
          coordinates[local] = T::variable( value, local );
          local_size = std::max( local_size, local + 1 );
        }
      }
    }
    return local_size;
  }

  template <typename AD>
  TRIBOL_FUTURE_HOST_DEVICE static void scatterGradient( const ContactDomainView& domain, ElementPair pair,
                                                         const AD& energy, int, Real* output )
  {
    const int elements[2] = { pair.mortar_element, pair.nonmortar_element };
    for ( int side = 0; side < 2; ++side ) {
      const int count = std::min( domain.numberOfElementDofs( elements[side] ), max_nodes_per_element );
      for ( int node = 0; node < count; ++node ) {
        const int dof = domain.elementDof( elements[side], node );
        for ( int component = 0; component < domain.dimension; ++component ) {
          const int local = ( side * max_nodes_per_element + node ) * 3 + component;
          AtomicAdd( output[component * domain.number_of_dofs + dof], detail::scalarValue( energy.derivative[local] ) );
        }
      }
    }
  }

  template <typename AD>
  TRIBOL_FUTURE_HOST_DEVICE static void scatterHessianVector( const ContactDomainView& domain, ElementPair pair,
                                                              const AD& energy, int, Real* output )
  {
    const int elements[2] = { pair.mortar_element, pair.nonmortar_element };
    for ( int side = 0; side < 2; ++side ) {
      const int count = std::min( domain.numberOfElementDofs( elements[side] ), max_nodes_per_element );
      for ( int node = 0; node < count; ++node ) {
        const int dof = domain.elementDof( elements[side], node );
        for ( int component = 0; component < domain.dimension; ++component ) {
          const int local = ( side * max_nodes_per_element + node ) * 3 + component;
          AtomicAdd( output[component * domain.number_of_dofs + dof], energy.derivative[local].derivative[0] );
        }
      }
    }
  }

  TRIBOL_FUTURE_HOST_DEVICE static void addVelocityResponse( const ContactDomainView& domain, ElementPair pair,
                                                             const Parameters& parameters, const Real* velocity,
                                                             const Real* velocity_direction, const Real* thickness,
                                                             const Real* modulus, Real* output )
  {
    if constexpr ( !capabilities.needs_velocity ) {
      return;
    }
    if ( velocity == nullptr && velocity_direction == nullptr ) {
      return;
    }
    const auto geometry = pairGeometry<Real>( domain, pair );
    const int elements[2] = { pair.mortar_element, pair.nonmortar_element };
    const int counts[2] = { std::min( domain.numberOfElementDofs( elements[0] ), max_nodes_per_element ),
                            std::min( domain.numberOfElementDofs( elements[1] ), max_nodes_per_element ) };
    std::array<Real, 3> relative_velocity{};
    const Real* field = velocity_direction != nullptr ? velocity_direction : velocity;
    for ( int component = 0; component < domain.dimension; ++component ) {
      for ( int node = 0; node < counts[0]; ++node ) {
        relative_velocity[component] +=
            field[component * domain.number_of_dofs + domain.elementDof( elements[0], node )] / counts[0];
      }
      for ( int node = 0; node < counts[1]; ++node ) {
        relative_velocity[component] -=
            field[component * domain.number_of_dofs + domain.elementDof( elements[1], node )] / counts[1];
      }
    }
    Real normal_rate{};
    for ( int component = 0; component < domain.dimension; ++component ) {
      normal_rate += relative_velocity[component] * geometry.normal[component];
    }
    const Real rate_coefficient =
        parameters.rate.rateCoefficient( normalStiffness( parameters.normal, thickness, modulus, pair ) );
    const Real normal_traction = rate_coefficient * std::min( normal_rate, 0.0 );
    std::array<Real, 3> traction{};
    for ( int component = 0; component < domain.dimension; ++component ) {
      traction[component] = normal_traction * geometry.normal[component];
    }
    parameters.tangential.addTraction( relative_velocity, normal_rate, geometry.normal, domain.dimension, traction );
    for ( int side = 0; side < 2; ++side ) {
      const Real sign = side == 0 ? 1.0 : -1.0;
      for ( int node = 0; node < counts[side]; ++node ) {
        const int dof = domain.elementDof( elements[side], node );
        for ( int component = 0; component < domain.dimension; ++component ) {
          AtomicAdd( output[component * domain.number_of_dofs + dof],
                     sign * geometry.measure * traction[component] / counts[side] );
        }
      }
    }
  }

  TRIBOL_FUTURE_HOST_DEVICE static void addVelocityCoordinateDerivative(
      const ContactDomainView& domain, ElementPair pair, const Parameters& parameters, const Real* velocity,
      const Real* thickness, const Real* modulus,
      const std::array<detail::Gradient<Real, 1>, max_local_coordinates>& coordinates, Real* output )
  {
    if constexpr ( !capabilities.needs_velocity ) {
      return;
    }
    if ( velocity == nullptr ) {
      return;
    }
    using Directional = detail::Gradient<Real, 1>;
    const auto geometry = pairGeometry( domain, pair, &coordinates );
    const int elements[2] = { pair.mortar_element, pair.nonmortar_element };
    const int counts[2] = { std::min( domain.numberOfElementDofs( elements[0] ), max_nodes_per_element ),
                            std::min( domain.numberOfElementDofs( elements[1] ), max_nodes_per_element ) };
    std::array<Real, 3> relative_velocity{};
    for ( int component = 0; component < domain.dimension; ++component ) {
      for ( int node = 0; node < counts[0]; ++node ) {
        relative_velocity[component] +=
            velocity[component * domain.number_of_dofs + domain.elementDof( elements[0], node )] / counts[0];
      }
      for ( int node = 0; node < counts[1]; ++node ) {
        relative_velocity[component] -=
            velocity[component * domain.number_of_dofs + domain.elementDof( elements[1], node )] / counts[1];
      }
    }
    Directional normal_rate{};
    for ( int component = 0; component < domain.dimension; ++component ) {
      normal_rate = normal_rate + Directional{ relative_velocity[component] } * geometry.normal[component];
    }
    const Real rate_coefficient =
        parameters.rate.rateCoefficient( normalStiffness( parameters.normal, thickness, modulus, pair ) );
    const Directional active_rate = detail::activeMinimum( normal_rate );
    std::array<Directional, 3> traction{};
    for ( int component = 0; component < domain.dimension; ++component ) {
      traction[component] = Directional{ rate_coefficient } * active_rate * geometry.normal[component];
    }
    std::array<Directional, 3> directional_velocity{};
    for ( int component = 0; component < domain.dimension; ++component ) {
      directional_velocity[component] = Directional{ relative_velocity[component] };
    }
    parameters.tangential.addTraction( directional_velocity, normal_rate, geometry.normal, domain.dimension, traction );
    for ( int side = 0; side < 2; ++side ) {
      const Real sign = side == 0 ? 1.0 : -1.0;
      for ( int node = 0; node < counts[side]; ++node ) {
        const int dof = domain.elementDof( elements[side], node );
        for ( int component = 0; component < domain.dimension; ++component ) {
          AtomicAdd( output[component * domain.number_of_dofs + dof],
                     sign * ( geometry.measure * traction[component] ).derivative[0] / counts[side] );
        }
      }
    }
  }

  void validateState( const ContactDomain& domain, const State& state ) const
  {
    if constexpr ( capabilities.needs_velocity ) {
      if ( state.velocity == nullptr || state.velocity->Size() != domain.coordinateSize() ) {
        throw std::invalid_argument( "CommonPlane requires a contact-space velocity vector." );
      }
    }
    if constexpr ( capabilities.needs_reference_field ) {
      if ( state.reference_coordinates == nullptr || state.reference_coordinates->Size() != domain.coordinateSize() ) {
        throw std::invalid_argument( "Tied CommonPlane response requires contact-space reference coordinates." );
      }
    }
    if constexpr ( capabilities.needs_material_fields ) {
      if ( state.element_thickness == nullptr || state.material_modulus == nullptr ||
           state.element_thickness->Size() != domain.numberOfElements() ||
           state.material_modulus->Size() != domain.numberOfElements() ) {
        throw std::invalid_argument( "MaterialKinematicPenalty requires element thickness and modulus fields." );
      }
    }
  }

  static void validateCoordinateVector( const ContactDomain& domain, const mfem::Vector& vector )
  {
    if ( vector.Size() != domain.coordinateSize() ) {
      throw std::invalid_argument( "CommonPlane coordinate vector has the wrong size." );
    }
  }

  static void setMemoryMode( Kinematics& kinematics, const InteractionStorage& interactions, bool use_device )
  {
    kinematics.gap.UseDevice( use_device );
    kinematics.pressure.UseDevice( use_device );
    kinematics.pressure_tangent.UseDevice( use_device );
    kinematics.potential_density.UseDevice( use_device );
    kinematics.pair_offsets.GetMemory().UseDevice( use_device );
    kinematics.quadrature_points.UseDevice( use_device );
    kinematics.quadrature_weights.UseDevice( use_device );
    kinematics.quadrature_gaps.UseDevice( use_device );
    kinematics.scalar_workspace.UseDevice( use_device );
    const_cast<mfem::Array<ElementPair>&>( interactions.pairs ).GetMemory().UseDevice( use_device );
    const_cast<InteractionColoring&>( interactions.coloring ).offsets.GetMemory().UseDevice( use_device );
    const_cast<InteractionColoring&>( interactions.coloring ).pair_indices.GetMemory().UseDevice( use_device );
  }

  static constexpr bool tiedNormal = NormalLaw::ties_normal;

  Parameters parameters_;
};

}  // namespace tribol::future

#endif
