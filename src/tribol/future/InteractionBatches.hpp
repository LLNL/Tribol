#ifndef TRIBOL_FUTURE_INTERACTIONBATCHES_HPP_
#define TRIBOL_FUTURE_INTERACTIONBATCHES_HPP_

#include "tribol/future/ContactDomain.hpp"
#include "tribol/future/Execution.hpp"

#include <algorithm>
#include <tuple>
#include <unordered_set>
#include <utility>
#include <vector>

namespace tribol::future {

struct InteractionBatch {
  Index begin;
  Index end;
  ElementTopology mortar_topology;
  ElementTopology nonmortar_topology;
  int mortar_order;
  int nonmortar_order;
};

struct InteractionColoring {
  mfem::Array<Index> offsets;
  mfem::Array<Index> pair_indices;
};

inline void buildInteractionBatches( const ContactDomain& domain, std::vector<ElementPair> pairs,
                                     mfem::Array<ElementPair>& ordered_pairs, mfem::Array<InteractionBatch>& batches )
{
  const auto view = domain.view( false );
  const auto key = [view]( ElementPair pair ) {
    return std::tuple{ view.topology( pair.mortar_element ),
                       view.topology( pair.nonmortar_element ),
                       view.orders[pair.mortar_element],
                       view.orders[pair.nonmortar_element],
                       pair.mortar_element,
                       pair.nonmortar_element };
  };
  std::sort( pairs.begin(), pairs.end(),
             [&]( ElementPair left, ElementPair right ) { return key( left ) < key( right ); } );

  ordered_pairs.SetSize( static_cast<int>( pairs.size() ) );
  std::vector<InteractionBatch> host_batches;
  for ( int pair_index = 0; pair_index < ordered_pairs.Size(); ++pair_index ) {
    const auto pair = pairs[static_cast<std::size_t>( pair_index )];
    ordered_pairs[pair_index] = pair;
    const InteractionBatch batch_key{ pair_index,
                                      pair_index + 1,
                                      view.topology( pair.mortar_element ),
                                      view.topology( pair.nonmortar_element ),
                                      view.orders[pair.mortar_element],
                                      view.orders[pair.nonmortar_element] };
    if ( host_batches.empty() || std::tie( host_batches.back().mortar_topology, host_batches.back().nonmortar_topology,
                                           host_batches.back().mortar_order, host_batches.back().nonmortar_order ) !=
                                     std::tie( batch_key.mortar_topology, batch_key.nonmortar_topology,
                                               batch_key.mortar_order, batch_key.nonmortar_order ) ) {
      host_batches.push_back( batch_key );
    } else {
      host_batches.back().end = pair_index + 1;
    }
  }

  batches.SetSize( static_cast<int>( host_batches.size() ) );
  for ( int batch = 0; batch < batches.Size(); ++batch ) {
    batches[batch] = host_batches[static_cast<std::size_t>( batch )];
  }
}

inline void buildInteractionColoring( const ContactDomain& domain, const mfem::Array<ElementPair>& pairs,
                                      InteractionColoring& coloring )
{
  std::vector<std::vector<Index>> colors;
  std::vector<std::unordered_set<Index>> color_dofs;
  const auto& element_offsets = domain.elementOffsets();
  const auto& element_dofs = domain.elementDofs();

  for ( Index pair_index = 0; pair_index < pairs.Size(); ++pair_index ) {
    const auto pair = pairs[pair_index];
    std::vector<Index> pair_dofs;
    for ( const Index element : { pair.mortar_element, pair.nonmortar_element } ) {
      for ( Index offset = element_offsets[element]; offset < element_offsets[element + 1]; ++offset ) {
        const Index dof = element_dofs[offset];
        if ( std::find( pair_dofs.begin(), pair_dofs.end(), dof ) == pair_dofs.end() ) {
          pair_dofs.push_back( dof );
        }
      }
    }

    std::size_t color = 0;
    for ( ; color < colors.size(); ++color ) {
      const bool conflict = std::any_of( pair_dofs.begin(), pair_dofs.end(),
                                         [&]( Index dof ) { return color_dofs[color].contains( dof ); } );
      if ( !conflict ) {
        break;
      }
    }
    if ( color == colors.size() ) {
      colors.emplace_back();
      color_dofs.emplace_back();
    }
    colors[color].push_back( pair_index );
    color_dofs[color].insert( pair_dofs.begin(), pair_dofs.end() );
  }

  coloring.offsets.SetSize( static_cast<int>( colors.size() ) + 1 );
  coloring.pair_indices.SetSize( pairs.Size() );
  Index offset = 0;
  coloring.offsets[0] = offset;
  for ( std::size_t color = 0; color < colors.size(); ++color ) {
    for ( const Index pair_index : colors[color] ) {
      coloring.pair_indices[offset++] = pair_index;
    }
    coloring.offsets[static_cast<int>( color ) + 1] = offset;
  }
}

template <ExecutionPolicy Execution, typename Body>
void forEachInteraction( const mfem::Array<InteractionBatch>& batches, const InteractionColoring& coloring, Body body )
{
  if constexpr ( Execution::deterministic ) {
    const Index* pair_indices = coloring.pair_indices.Read( Execution::uses_device );
    using Backend = typename Execution::backend_type;
    for ( int color = 0; color + 1 < coloring.offsets.Size(); ++color ) {
      const Index begin = coloring.offsets[color];
      const Index end = coloring.offsets[color + 1];
      Backend::forAll( end - begin, [=] MFEM_HOST_DEVICE( int index ) { body( pair_indices[begin + index] ); } );
    }
  } else {
    for ( int batch_index = 0; batch_index < batches.Size(); ++batch_index ) {
      const Index begin = batches[batch_index].begin;
      const Index end = batches[batch_index].end;
      Execution::forAll( end - begin, [=] MFEM_HOST_DEVICE( int index ) { body( begin + index ); } );
    }
  }
}

template <ExecutionPolicy Execution>
TRIBOL_FUTURE_HOST_DEVICE void storeInteractionScalar( Real* workspace, int pair_index, Real value )
{
  if constexpr ( Execution::deterministic ) {
    workspace[pair_index + 1] = value;
  } else {
    AtomicAdd( workspace[0], value );
  }
}

template <ExecutionPolicy Execution>
void finalizeInteractionScalar( mfem::Vector& workspace, int pair_count )
{
  if constexpr ( Execution::deterministic ) {
    Real* values = workspace.ReadWrite( Execution::uses_device );
    using Backend = typename Execution::backend_type;
    Backend::forAll( 1, [=] MFEM_HOST_DEVICE( int ) {
      Real sum{};
      for ( int pair_index = 0; pair_index < pair_count; ++pair_index ) {
        sum += values[pair_index + 1];
      }
      values[0] = sum;
    } );
  }
}

}  // namespace tribol::future

#endif
