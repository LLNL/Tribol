#ifndef TRIBOL_FUTURE_SEARCH_HPP_
#define TRIBOL_FUTURE_SEARCH_HPP_

#include "tribol/future/ContactDomain.hpp"
#include "tribol/future/Execution.hpp"

#include "axom/core/execution/execution_space.hpp"
#include "axom/core/execution/for_all.hpp"
#include "axom/core/memory_management.hpp"
#include "axom/primal.hpp"
#include "axom/spin/BVH.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

namespace tribol::future {

namespace detail {

template <int Dimension>
axom::primal::BoundingBox<Real, Dimension> elementBoundingBox( const ContactDomainView& domain, int element,
                                                               Real expansion )
{
  axom::primal::BoundingBox<Real, Dimension> box;
  for ( int local_dof = 0; local_dof < domain.numberOfElementDofs( element ); ++local_dof ) {
    const int dof = domain.elementDof( element, local_dof );
    axom::primal::Point<Real, Dimension> point;
    for ( int component = 0; component < Dimension; ++component ) {
      point[component] = domain.coordinate( component, dof );
    }
    box.addPoint( point );
  }
  if ( box.isValid() && expansion > 0.0 ) {
    auto lower = box.getMin();
    auto upper = box.getMax();
    for ( int component = 0; component < Dimension; ++component ) {
      lower[component] -= expansion;
      upper[component] += expansion;
    }
    box = axom::primal::BoundingBox<Real, Dimension>( lower, upper );
  }
  return box;
}

TRIBOL_FUTURE_HOST_DEVICE inline bool shareDof( const ContactDomainView& domain, int first, int second )
{
  for ( int i = 0; i < domain.numberOfElementDofs( first ); ++i ) {
    for ( int j = 0; j < domain.numberOfElementDofs( second ); ++j ) {
      if ( domain.elementDof( first, i ) == domain.elementDof( second, j ) ) {
        return true;
      }
    }
  }
  return false;
}

inline bool withinAdjacencyDepth( const ContactDomainView& domain, int first, int second, int adjacency_depth )
{
  if ( adjacency_depth <= 0 ) {
    return false;
  }
  std::vector<unsigned char> visited( static_cast<std::size_t>( domain.number_of_elements ), 0 );
  std::vector<int> frontier{ first };
  visited[static_cast<std::size_t>( first )] = 1;
  for ( int depth = 0; depth < adjacency_depth; ++depth ) {
    std::vector<int> next;
    for ( const int element : frontier ) {
      for ( int neighbor = 0; neighbor < domain.number_of_elements; ++neighbor ) {
        if ( visited[static_cast<std::size_t>( neighbor )] != 0 || !shareDof( domain, element, neighbor ) ) {
          continue;
        }
        if ( neighbor == second ) {
          return true;
        }
        visited[static_cast<std::size_t>( neighbor )] = 1;
        next.push_back( neighbor );
      }
    }
    frontier = std::move( next );
    if ( frontier.empty() ) {
      break;
    }
  }
  return false;
}

TRIBOL_FUTURE_HOST_DEVICE inline ElementPair canonicalPair( const ContactDomainView& domain, int first, int second )
{
  const bool first_precedes = domain.owner_ranks[first] < domain.owner_ranks[second] ||
                              ( domain.owner_ranks[first] == domain.owner_ranks[second] &&
                                domain.owner_indices[first] < domain.owner_indices[second] );
  return first_precedes ? ElementPair{ first, second } : ElementPair{ second, first };
}

template <int Dimension>
void cartesianCandidates( const ContactDomainView& domain, Real expansion, bool self_contact, int adjacency_depth,
                          std::vector<ElementPair>& result )
{
  for ( int first = 0; first < domain.number_of_elements; ++first ) {
    if ( !self_contact && domain.side( first ) != SurfaceSide::Mortar ) {
      continue;
    }
    const auto first_box = elementBoundingBox<Dimension>( domain, first, expansion );
    for ( int second = self_contact ? first + 1 : 0; second < domain.number_of_elements; ++second ) {
      if ( !self_contact && domain.side( second ) != SurfaceSide::Nonmortar ) {
        continue;
      }
      if ( first == second ) {
        continue;
      }
      if ( self_contact && withinAdjacencyDepth( domain, first, second, adjacency_depth ) ) {
        continue;
      }
      if ( first_box.intersectsWith( elementBoundingBox<Dimension>( domain, second, expansion ) ) ) {
        const auto pair = self_contact ? canonicalPair( domain, first, second ) : ElementPair{ first, second };
        if ( domain.ghost_elements[pair.mortar_element] == 0 ) {
          result.push_back( pair );
        }
      }
    }
  }
}

template <int Dimension>
void bvhCandidates( const ContactDomainView& domain, Real expansion, bool self_contact, int adjacency_depth,
                    std::vector<ElementPair>& result )
{
  using Box = axom::primal::BoundingBox<Real, Dimension>;
  std::vector<int> indexed_elements;
  std::vector<int> query_elements;
  std::vector<Box> indexed_boxes;
  std::vector<Box> query_boxes;

  for ( int element = 0; element < domain.number_of_elements; ++element ) {
    const auto box = elementBoundingBox<Dimension>( domain, element, expansion );
    if ( self_contact || domain.side( element ) == SurfaceSide::Nonmortar ) {
      indexed_elements.push_back( element );
      indexed_boxes.push_back( box );
    }
    if ( self_contact || domain.side( element ) == SurfaceSide::Mortar ) {
      query_elements.push_back( element );
      query_boxes.push_back( box );
    }
  }

  if ( indexed_boxes.empty() || query_boxes.empty() ) {
    return;
  }

  axom::spin::BVH<Dimension, axom::SEQ_EXEC, Real> bvh;
  bvh.initialize( indexed_boxes.data(), static_cast<axom::IndexType>( indexed_boxes.size() ) );
  axom::Array<axom::IndexType> offsets( query_boxes.size(), query_boxes.size() );
  axom::Array<axom::IndexType> counts( query_boxes.size(), query_boxes.size() );
  axom::Array<axom::IndexType> candidates;
  bvh.findBoundingBoxes( offsets.view(), counts.view(), candidates, query_boxes.size(), query_boxes.data() );

  for ( int query = 0; query < static_cast<int>( query_elements.size() ); ++query ) {
    const int first = query_elements[query];
    for ( int candidate = 0; candidate < counts[query]; ++candidate ) {
      const int second = indexed_elements[candidates[offsets[query] + candidate]];
      if ( first == second ) {
        continue;
      }
      if ( self_contact ) {
        const auto pair = canonicalPair( domain, first, second );
        if ( withinAdjacencyDepth( domain, pair.mortar_element, pair.nonmortar_element, adjacency_depth ) ) {
          continue;
        }
        if ( domain.ghost_elements[pair.mortar_element] == 0 ) {
          result.push_back( pair );
        }
      } else {
        if ( domain.ghost_elements[first] == 0 ) {
          result.push_back( { first, second } );
        }
      }
    }
  }

  std::sort( result.begin(), result.end(), []( const ElementPair& left, const ElementPair& right ) {
    return std::pair{ left.mortar_element, left.nonmortar_element } <
           std::pair{ right.mortar_element, right.nonmortar_element };
  } );
  result.erase( std::unique( result.begin(), result.end() ), result.end() );
}

template <int Dimension, typename AxomExecution>
void deviceBvhCandidates( const ContactDomainView& host_domain, Real expansion, bool self_contact, int adjacency_depth,
                          std::vector<ElementPair>& result )
{
  using Box = axom::primal::BoundingBox<Real, Dimension>;
  using IndexArray = axom::Array<axom::IndexType>;
  const int host_allocator = axom::execution_space<axom::SEQ_EXEC>::allocatorID();
  const int device_allocator = axom::execution_space<AxomExecution>::allocatorID();

  std::vector<int> indexed_elements;
  std::vector<int> query_elements;
  axom::Array<Box> indexed_host;
  axom::Array<Box> query_host;
  for ( int element = 0; element < host_domain.number_of_elements; ++element ) {
    const auto box = elementBoundingBox<Dimension>( host_domain, element, expansion );
    if ( self_contact || host_domain.side( element ) == SurfaceSide::Nonmortar ) {
      indexed_elements.push_back( element );
      indexed_host.push_back( box );
    }
    if ( self_contact || host_domain.side( element ) == SurfaceSide::Mortar ) {
      query_elements.push_back( element );
      query_host.push_back( box );
    }
  }
  if ( indexed_host.empty() || query_host.empty() ) {
    return;
  }

  axom::Array<Box> indexed_device( indexed_host, device_allocator );
  axom::Array<Box> query_device( query_host, device_allocator );
  axom::spin::BVH<Dimension, AxomExecution, Real> bvh;
  bvh.setAllocatorID( device_allocator );
  bvh.initialize( indexed_device.view(), indexed_device.size() );
  IndexArray offsets_device( query_device.size(), query_device.size(), device_allocator );
  IndexArray counts_device( query_device.size(), query_device.size(), device_allocator );
  IndexArray candidates_device( 0, 0, device_allocator );
  bvh.findBoundingBoxes( offsets_device.view(), counts_device.view(), candidates_device, query_device.size(),
                         query_device.view() );

  IndexArray offsets_host( offsets_device, host_allocator );
  IndexArray counts_host( counts_device, host_allocator );
  IndexArray candidates_host( candidates_device, host_allocator );
  for ( int query = 0; query < static_cast<int>( query_elements.size() ); ++query ) {
    const int first = query_elements[static_cast<std::size_t>( query )];
    for ( int candidate = 0; candidate < counts_host[query]; ++candidate ) {
      const int second = indexed_elements[static_cast<std::size_t>( candidates_host[offsets_host[query] + candidate] )];
      if ( first == second ) {
        continue;
      }
      if ( self_contact ) {
        const auto pair = canonicalPair( host_domain, first, second );
        if ( withinAdjacencyDepth( host_domain, pair.mortar_element, pair.nonmortar_element, adjacency_depth ) ) {
          continue;
        }
        if ( host_domain.ghost_elements[pair.mortar_element] == 0 ) {
          result.push_back( pair );
        }
      } else {
        if ( host_domain.ghost_elements[first] == 0 ) {
          result.push_back( { first, second } );
        }
      }
    }
  }
  std::sort( result.begin(), result.end(), []( const ElementPair& left, const ElementPair& right ) {
    return std::pair{ left.mortar_element, left.nonmortar_element } <
           std::pair{ right.mortar_element, right.nonmortar_element };
  } );
  result.erase( std::unique( result.begin(), result.end() ), result.end() );
}

inline void assignCandidates( const std::vector<ElementPair>& source, CandidatePairs& destination )
{
  destination.pairs.SetSize( static_cast<int>( source.size() ) );
  for ( int i = 0; i < destination.pairs.Size(); ++i ) {
    destination.pairs[i] = source[static_cast<std::size_t>( i )];
  }
}

}  // namespace detail

class CartesianProductSearch {
 public:
  struct Parameters {
    Real expansion{ 0.0 };
  };

  CartesianProductSearch() = default;
  explicit CartesianProductSearch( Parameters parameters ) : parameters_( parameters ) {}

  void findCandidates( const ContactDomain& domain, CandidatePairs& candidates, bool self_contact,
                       int adjacency_depth ) const
  {
    const auto view = domain.view( false );
    std::vector<ElementPair> result;
    if ( view.dimension == 2 ) {
      detail::cartesianCandidates<2>( view, parameters_.expansion, self_contact, adjacency_depth, result );
    } else {
      detail::cartesianCandidates<3>( view, parameters_.expansion, self_contact, adjacency_depth, result );
    }
    detail::assignCandidates( result, candidates );
  }

 private:
  Parameters parameters_;
};

class BvhSearch {
 public:
  struct Parameters {
    Real expansion{ 0.0 };
  };

  BvhSearch() = default;
  explicit BvhSearch( Parameters parameters ) : parameters_( parameters ) {}

  void findCandidates( const ContactDomain& domain, CandidatePairs& candidates, bool self_contact,
                       int adjacency_depth ) const
  {
    const auto view = domain.view( false );
    std::vector<ElementPair> result;
    if ( view.dimension == 2 ) {
      detail::bvhCandidates<2>( view, parameters_.expansion, self_contact, adjacency_depth, result );
    } else {
      detail::bvhCandidates<3>( view, parameters_.expansion, self_contact, adjacency_depth, result );
    }
    detail::assignCandidates( result, candidates );
  }

  template <ExecutionPolicy Execution>
  void findCandidates( const ContactDomain& domain, CandidatePairs& candidates, bool self_contact,
                       int adjacency_depth ) const
  {
    if constexpr ( Execution::uses_hip ) {
#if defined( AXOM_USE_HIP ) && defined( MFEM_USE_HIP )
      const auto host_view = domain.view( false );
      std::vector<ElementPair> result;
      if ( host_view.dimension == 2 ) {
        detail::deviceBvhCandidates<2, axom::HIP_EXEC<256>>( host_view, parameters_.expansion, self_contact,
                                                             adjacency_depth, result );
      } else {
        detail::deviceBvhCandidates<3, axom::HIP_EXEC<256>>( host_view, parameters_.expansion, self_contact,
                                                             adjacency_depth, result );
      }
      detail::assignCandidates( result, candidates );
#else
      static_assert( !Execution::uses_device, "HIP BVH search requires HIP-enabled Axom and MFEM." );
#endif
    } else if constexpr ( Execution::uses_cuda ) {
#if defined( AXOM_USE_CUDA ) && defined( MFEM_USE_CUDA )
      const auto host_view = domain.view( false );
      std::vector<ElementPair> result;
      if ( host_view.dimension == 2 ) {
        detail::deviceBvhCandidates<2, axom::CUDA_EXEC<256>>( host_view, parameters_.expansion, self_contact,
                                                              adjacency_depth, result );
      } else {
        detail::deviceBvhCandidates<3, axom::CUDA_EXEC<256>>( host_view, parameters_.expansion, self_contact,
                                                              adjacency_depth, result );
      }
      detail::assignCandidates( result, candidates );
#else
      static_assert( !Execution::uses_device, "CUDA BVH search requires CUDA-enabled Axom and MFEM." );
#endif
    } else {
      findCandidates( domain, candidates, self_contact, adjacency_depth );
    }
  }

 private:
  Parameters parameters_;
};

}  // namespace tribol::future

#endif
