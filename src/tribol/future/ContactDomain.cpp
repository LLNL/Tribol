#include "tribol/future/ContactDomain.hpp"

#include "redecomp/RedecompMesh.hpp"

#include <algorithm>
#include <stdexcept>
#include <vector>

namespace tribol::future {

namespace {

int decodeDof( int dof ) { return dof >= 0 ? dof : -1 - dof; }

}  // namespace

void ContactDomain::rebuild( const mfem::Mesh& mesh, const mfem::GridFunction& coordinates,
                             const std::vector<int>& mortar_attributes, const std::vector<int>& nonmortar_attributes,
                             bool self_contact )
{
  const auto* coordinate_space = coordinates.FESpace();
  if ( coordinate_space == nullptr || coordinate_space->GetMesh() != &mesh ) {
    throw std::invalid_argument( "ContactDomain coordinates must be defined on its mesh." );
  }

  dimension_ = coordinate_space->GetVDim();
  number_of_dofs_ = coordinate_space->GetNDofs();
  if ( dimension_ != 2 && dimension_ != 3 ) {
    throw std::invalid_argument( "ContactDomain supports only two- and three-dimensional coordinates." );
  }

  const int number_of_elements = mesh.GetNE();
  element_offsets_.SetSize( number_of_elements + 1 );
  element_attributes_.SetSize( number_of_elements );
  element_sides_.SetSize( number_of_elements );
  element_topologies_.SetSize( number_of_elements );
  element_orders_.SetSize( number_of_elements );
  owner_ranks_.SetSize( number_of_elements );
  owner_indices_.SetSize( number_of_elements );
  ghost_elements_.SetSize( number_of_elements );

  std::vector<Index> dofs;
  std::vector<Real> reference_coordinates;
  dofs.reserve( static_cast<std::size_t>( number_of_elements ) * 4 );
  reference_coordinates.reserve( static_cast<std::size_t>( number_of_elements ) * 8 );
  element_offsets_[0] = 0;

  const auto* redecomp_mesh = dynamic_cast<const redecomp::RedecompMesh*>( &mesh );
  const auto* owner_offsets = redecomp_mesh ? &redecomp_mesh->getRedecompToParentElemOffsets() : nullptr;
  const auto* ghost_elements = redecomp_mesh ? &redecomp_mesh->getRedecompToParentGhostElems() : nullptr;
  int owner_rank = 0;

  for ( int element = 0; element < number_of_elements; ++element ) {
    mfem::Array<int> element_dofs;
    coordinate_space->GetElementDofs( element, element_dofs );
    const auto& nodes = coordinate_space->GetFE( element )->GetNodes();
    if ( nodes.GetNPoints() != element_dofs.Size() ) {
      throw std::runtime_error( "ContactDomain requires a nodal coordinate finite-element space." );
    }
    for ( int local_dof = 0; local_dof < element_dofs.Size(); ++local_dof ) {
      dofs.push_back( static_cast<Index>( decodeDof( element_dofs[local_dof] ) ) );
      const auto& point = nodes.IntPoint( local_dof );
      reference_coordinates.push_back( point.x );
      reference_coordinates.push_back( point.y );
    }
    element_offsets_[element + 1] = static_cast<Index>( dofs.size() );

    const int attribute = mesh.GetAttribute( element );
    element_attributes_[element] = attribute;
    if ( self_contact || contains( mortar_attributes, attribute ) ) {
      element_sides_[element] = static_cast<unsigned char>( SurfaceSide::Mortar );
    } else if ( contains( nonmortar_attributes, attribute ) ) {
      element_sides_[element] = static_cast<unsigned char>( SurfaceSide::Nonmortar );
    } else {
      throw std::runtime_error( "A contact-domain element has an unclassified boundary attribute." );
    }
    element_topologies_[element] = static_cast<unsigned char>( topology( mesh.GetElementBaseGeometry( element ) ) );
    element_orders_[element] = coordinate_space->GetFE( element )->GetOrder();

    if ( owner_offsets != nullptr ) {
      while ( owner_rank + 1 < owner_offsets->size() && element >= ( *owner_offsets )[owner_rank + 1] ) {
        ++owner_rank;
      }
      owner_ranks_[element] = owner_rank;
      owner_indices_[element] = element - ( *owner_offsets )[owner_rank];
    } else {
      owner_ranks_[element] = 0;
      owner_indices_[element] = element;
    }
    ghost_elements_[element] = 0;
    if ( ghost_elements != nullptr ) {
      const auto& owner_ghosts = ( *ghost_elements )[owner_rank];
      ghost_elements_[element] = std::find( owner_ghosts.begin(), owner_ghosts.end(), element ) != owner_ghosts.end();
    }
  }

  element_dofs_.SetSize( static_cast<int>( dofs.size() ) );
  for ( int i = 0; i < element_dofs_.Size(); ++i ) {
    element_dofs_[i] = dofs[static_cast<std::size_t>( i )];
  }
  reference_coordinates_.SetSize( static_cast<int>( reference_coordinates.size() ) );
  for ( int i = 0; i < reference_coordinates_.Size(); ++i ) {
    reference_coordinates_[i] = reference_coordinates[static_cast<std::size_t>( i )];
  }

  coordinates_.SetSize( dimension_ * number_of_dofs_ );
  coordinates_.UseDevice( coordinates.UseDevice() );
  updateCoordinates( coordinates );

  const bool use_device = coordinates.UseDevice();
  element_offsets_.GetMemory().UseDevice( use_device );
  element_dofs_.GetMemory().UseDevice( use_device );
  element_attributes_.GetMemory().UseDevice( use_device );
  element_sides_.GetMemory().UseDevice( use_device );
  element_topologies_.GetMemory().UseDevice( use_device );
  element_orders_.GetMemory().UseDevice( use_device );
  reference_coordinates_.UseDevice( use_device );
  owner_ranks_.GetMemory().UseDevice( use_device );
  owner_indices_.GetMemory().UseDevice( use_device );
  ghost_elements_.GetMemory().UseDevice( use_device );
}

void ContactDomain::updateCoordinates( const mfem::GridFunction& coordinates )
{
  const auto* coordinate_space = coordinates.FESpace();
  if ( coordinate_space == nullptr || coordinate_space->GetNDofs() != number_of_dofs_ ||
       coordinate_space->GetVDim() != dimension_ ) {
    throw std::invalid_argument( "Updated contact coordinates do not match the existing contact domain." );
  }

  const bool use_device = coordinates.UseDevice();
  const Real* source = coordinates.Read( use_device );
  Real* destination = coordinates_.Write( use_device );
  const int number_of_dofs = number_of_dofs_;
  const int dimension = dimension_;
  const auto ordering = coordinate_space->GetOrdering();
  mfem::forall_switch( use_device, dimension_ * number_of_dofs_, [=] MFEM_HOST_DEVICE( int index ) {
    const int component = index / number_of_dofs;
    const int dof = index - component * number_of_dofs;
    const int source_index =
        ordering == mfem::Ordering::byNODES ? component * number_of_dofs + dof : dof * dimension + component;
    destination[index] = source[source_index];
  } );
}

ContactDomainView ContactDomain::view( bool on_device ) const
{
  return { coordinates_.Read( on_device ),
           element_offsets_.Read( on_device ),
           element_dofs_.Read( on_device ),
           element_attributes_.Read( on_device ),
           element_sides_.Read( on_device ),
           element_topologies_.Read( on_device ),
           element_orders_.Read( on_device ),
           reference_coordinates_.Read( on_device ),
           owner_ranks_.Read( on_device ),
           owner_indices_.Read( on_device ),
           ghost_elements_.Read( on_device ),
           dimension_,
           number_of_dofs_,
           numberOfElements() };
}

ElementTopology ContactDomain::topology( mfem::Geometry::Type geometry )
{
  switch ( geometry ) {
    case mfem::Geometry::SEGMENT:
      return ElementTopology::Segment;
    case mfem::Geometry::TRIANGLE:
      return ElementTopology::Triangle;
    case mfem::Geometry::SQUARE:
      return ElementTopology::Quadrilateral;
    default:
      return ElementTopology::Unsupported;
  }
}

bool ContactDomain::contains( const std::vector<int>& values, int value )
{
  return std::find( values.begin(), values.end(), value ) != values.end();
}

}  // namespace tribol::future
