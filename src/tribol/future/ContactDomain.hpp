#ifndef TRIBOL_FUTURE_CONTACTDOMAIN_HPP_
#define TRIBOL_FUTURE_CONTACTDOMAIN_HPP_

#include "tribol/future/Types.hpp"

#include "mfem.hpp"

#include <type_traits>

namespace tribol::future {

struct ContactDomainView {
  const Real* coordinates{};
  const Index* element_offsets{};
  const Index* element_dofs{};
  const int* attributes{};
  const unsigned char* sides{};
  const unsigned char* topologies{};
  const int* orders{};
  const Real* reference_coordinates{};
  const int* owner_ranks{};
  const int* owner_indices{};
  const unsigned char* ghost_elements{};
  int dimension{};
  int number_of_dofs{};
  int number_of_elements{};

  TRIBOL_FUTURE_HOST_DEVICE Real coordinate( int component, int dof ) const
  {
    return coordinates[component * number_of_dofs + dof];
  }

  TRIBOL_FUTURE_HOST_DEVICE int numberOfElementDofs( int element ) const
  {
    return element_offsets[element + 1] - element_offsets[element];
  }

  TRIBOL_FUTURE_HOST_DEVICE int elementDof( int element, int local_dof ) const
  {
    return element_dofs[element_offsets[element] + local_dof];
  }

  TRIBOL_FUTURE_HOST_DEVICE Real referenceCoordinate( int element, int local_dof, int component = 0 ) const
  {
    return reference_coordinates[2 * ( element_offsets[element] + local_dof ) + component];
  }

  TRIBOL_FUTURE_HOST_DEVICE SurfaceSide side( int element ) const { return static_cast<SurfaceSide>( sides[element] ); }

  TRIBOL_FUTURE_HOST_DEVICE ElementTopology topology( int element ) const
  {
    return static_cast<ElementTopology>( topologies[element] );
  }
};

static_assert( std::is_trivially_copyable_v<ContactDomainView> );

class ContactDomain {
 public:
  ContactDomain() = default;

  void rebuild( const mfem::Mesh& mesh, const mfem::GridFunction& coordinates,
                const std::vector<int>& mortar_attributes, const std::vector<int>& nonmortar_attributes,
                bool self_contact );

  void updateCoordinates( const mfem::GridFunction& coordinates );

  ContactDomainView view( bool on_device ) const;

  int dimension() const { return dimension_; }
  int numberOfDofs() const { return number_of_dofs_; }
  int numberOfElements() const { return element_attributes_.Size(); }
  int coordinateSize() const { return coordinates_.Size(); }

  const mfem::Vector& coordinates() const { return coordinates_; }
  const mfem::Array<Index>& elementOffsets() const { return element_offsets_; }
  const mfem::Array<Index>& elementDofs() const { return element_dofs_; }
  const mfem::Array<int>& elementAttributes() const { return element_attributes_; }
  const mfem::Array<unsigned char>& elementSides() const { return element_sides_; }
  const mfem::Array<unsigned char>& elementTopologies() const { return element_topologies_; }
  const mfem::Array<int>& elementOrders() const { return element_orders_; }
  const mfem::Vector& referenceCoordinates() const { return reference_coordinates_; }

 private:
  static ElementTopology topology( mfem::Geometry::Type geometry );
  static bool contains( const std::vector<int>& values, int value );

  int dimension_{};
  int number_of_dofs_{};
  mfem::Vector coordinates_;
  mfem::Array<Index> element_offsets_;
  mfem::Array<Index> element_dofs_;
  mfem::Array<int> element_attributes_;
  mfem::Array<unsigned char> element_sides_;
  mfem::Array<unsigned char> element_topologies_;
  mfem::Array<int> element_orders_;
  mfem::Vector reference_coordinates_;
  mfem::Array<int> owner_ranks_;
  mfem::Array<int> owner_indices_;
  mfem::Array<unsigned char> ghost_elements_;
};

struct CandidatePairs {
  mfem::Array<ElementPair> pairs;

  int size() const { return pairs.Size(); }
  bool empty() const { return pairs.Size() == 0; }
};

}  // namespace tribol::future

#endif
