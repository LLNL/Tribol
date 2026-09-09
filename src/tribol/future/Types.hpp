#ifndef TRIBOL_FUTURE_TYPES_HPP_
#define TRIBOL_FUTURE_TYPES_HPP_

#include "tribol/future/Config.hpp"

#include <initializer_list>
#include <stdexcept>
#include <utility>
#include <vector>

namespace tribol::future {

enum class SurfaceSide : unsigned char
{
  Mortar = 0,
  Nonmortar = 1
};

enum class ElementTopology : unsigned char
{
  Segment,
  Triangle,
  Quadrilateral,
  Unsupported
};

struct MethodCapabilities {
  bool needs_dual_space{};
  bool needs_velocity{};
  bool needs_reference_field{};
  bool needs_material_fields{};
  bool has_history{};
  bool has_energy{};
  bool supports_self_contact{};
  bool supports_assembled_jacobian{};
};

struct PairedSurfaces {
  std::vector<int> mortar_attributes;
  std::vector<int> nonmortar_attributes;

  PairedSurfaces( std::initializer_list<int> mortar, std::initializer_list<int> nonmortar )
      : mortar_attributes( mortar ), nonmortar_attributes( nonmortar )
  {
    if ( mortar_attributes.empty() || nonmortar_attributes.empty() ) {
      throw std::invalid_argument( "PairedSurfaces requires two nonempty boundary-attribute sets." );
    }
  }

  PairedSurfaces( std::vector<int> mortar, std::vector<int> nonmortar )
      : mortar_attributes( std::move( mortar ) ), nonmortar_attributes( std::move( nonmortar ) )
  {
    if ( mortar_attributes.empty() || nonmortar_attributes.empty() ) {
      throw std::invalid_argument( "PairedSurfaces requires two nonempty boundary-attribute sets." );
    }
  }
};

struct SelfContactSurface {
  std::vector<int> attributes;
  int adjacency_exclusion_depth{ 1 };

  SelfContactSurface( std::initializer_list<int> attributes_, int exclusion_depth = 1 )
      : attributes( attributes_ ), adjacency_exclusion_depth( exclusion_depth )
  {
    validate();
  }

  explicit SelfContactSurface( std::vector<int> attributes_, int exclusion_depth = 1 )
      : attributes( std::move( attributes_ ) ), adjacency_exclusion_depth( exclusion_depth )
  {
    validate();
  }

 private:
  void validate() const
  {
    if ( attributes.empty() ) {
      throw std::invalid_argument( "SelfContactSurface requires a nonempty boundary-attribute set." );
    }
    if ( adjacency_exclusion_depth < 0 ) {
      throw std::invalid_argument( "Self-contact adjacency exclusion depth cannot be negative." );
    }
  }
};

struct ElementPair {
  Index mortar_element;
  Index nonmortar_element;
};

inline constexpr bool operator==( ElementPair left, ElementPair right )
{
  return left.mortar_element == right.mortar_element && left.nonmortar_element == right.nonmortar_element;
}

struct GeometryVersion {
  std::uint64_t value{};

  GeometryVersion& operator++()
  {
    ++value;
    return *this;
  }

  friend constexpr bool operator==( GeometryVersion, GeometryVersion ) = default;
};

struct InteractionVersion {
  std::uint64_t value{};

  InteractionVersion& operator++()
  {
    ++value;
    return *this;
  }

  friend constexpr bool operator==( InteractionVersion, InteractionVersion ) = default;
};

}  // namespace tribol::future

#endif
