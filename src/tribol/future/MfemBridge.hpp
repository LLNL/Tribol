#ifndef TRIBOL_FUTURE_MFEMBRIDGE_HPP_
#define TRIBOL_FUTURE_MFEMBRIDGE_HPP_

#include "tribol/future/ContactDomain.hpp"
#include "tribol/future/ContactRestrictionOperator.hpp"
#include "tribol/future/SurfaceRestrictionOperator.hpp"
#include "tribol/future/SurfaceDiscretizations.hpp"
#include "tribol/future/Types.hpp"

#include "mfem.hpp"

#include "redecomp/RedecompMesh.hpp"
#include "redecomp/RedecompTransfer.hpp"

#include <memory>
#include <variant>
#include <vector>

namespace tribol::future {

struct MfemBridgeParameters {
  Real ghost_distance{ -1.0 };
  int redecomp_ranks{};
  bool use_device{ false };
};

class MfemBridge {
 public:
  template <typename Surface>
  MfemBridge( std::type_identity<Surface>, const mfem::ParMesh& parent_mesh, const mfem::ParGridFunction& coordinates,
              const PairedSurfaces& surfaces, const typename Surface::Parameters& surface_parameters,
              MfemBridgeParameters parameters = {} )
      : parent_mesh_( parent_mesh ),
        coordinates_( &coordinates ),
        parent_mesh_sequence_( parent_mesh.GetSequence() ),
        parent_space_sequence_( coordinates.ParFESpace()->GetSequence() ),
        mortar_attributes_( surfaces.mortar_attributes ),
        nonmortar_attributes_( surfaces.nonmortar_attributes ),
        parameters_( parameters )
  {
    Surface::validate( coordinates, surface_parameters );
    configureSurface<Surface>( surface_parameters );
    rebuildGeometry( coordinates );
  }

  template <typename Surface>
  MfemBridge( std::type_identity<Surface>, const mfem::ParMesh& parent_mesh, const mfem::ParGridFunction& coordinates,
              const SelfContactSurface& surface, const typename Surface::Parameters& surface_parameters,
              MfemBridgeParameters parameters = {} )
      : parent_mesh_( parent_mesh ),
        coordinates_( &coordinates ),
        parent_mesh_sequence_( parent_mesh.GetSequence() ),
        parent_space_sequence_( coordinates.ParFESpace()->GetSequence() ),
        mortar_attributes_( surface.attributes ),
        self_contact_( true ),
        adjacency_exclusion_depth_( surface.adjacency_exclusion_depth ),
        parameters_( parameters )
  {
    Surface::validate( coordinates, surface_parameters );
    configureSurface<Surface>( surface_parameters );
    rebuildGeometry( coordinates );
  }

  void updateGeometry( const mfem::ParGridFunction& coordinates );
  void rebuildGeometry( const mfem::ParGridFunction& coordinates );

  const ContactDomain& domain() const { return domain_; }
  ContactDomain& domain() { return domain_; }
  const ContactRestrictionOperator& coordinateRestriction() const { return *coordinate_restriction_; }
  const SurfaceRestrictionOperator& scalarRestriction() const { return *scalar_restriction_; }
  const mfem::ParFiniteElementSpace& parentCoordinateSpace() const { return *coordinates_->ParFESpace(); }
  const mfem::ParFiniteElementSpace& surfaceCoordinateSpace() const { return *surface_coordinate_space_; }
  const mfem::ParFiniteElementSpace& surfaceScalarSpace() const { return *surface_scalar_space_; }
  const mfem::FiniteElementSpace& contactCoordinateSpace() const { return *redecomp_coordinates_->FESpace(); }
  const mfem::FiniteElementSpace& contactScalarSpace() const { return *redecomp_scalar_space_; }
  const mfem::GridFunction& contactCoordinates() const { return *redecomp_coordinates_; }

  GeometryVersion geometryVersion() const { return geometry_version_; }
  bool topologyIsCurrent() const
  {
    return parent_mesh_.GetSequence() == parent_mesh_sequence_ &&
           coordinates_->ParFESpace()->GetSequence() == parent_space_sequence_;
  }
  bool interactionsValid() const { return interactions_valid_; }
  void markInteractionsUpdated()
  {
    interactions_valid_ = true;
    ++interaction_version_;
  }
  void invalidateInteractions() { interactions_valid_ = false; }
  InteractionVersion interactionVersion() const { return interaction_version_; }

  bool selfContact() const { return self_contact_; }
  int adjacencyExclusionDepth() const { return adjacency_exclusion_depth_; }
  bool usesLowOrderRefinement() const { return use_lor_; }
  int maximumSurfaceOrder() const { return maximum_surface_order_; }

 private:
  template <typename Surface>
  void configureSurface( const typename Surface::Parameters& parameters )
  {
    use_lor_ = Surface::is_low_order_refined;
    if constexpr ( Surface::is_low_order_refined ) {
      lor_factor_ = parameters.refinement_factor;
    } else {
      maximum_order_ = Surface::max_order;
    }
  }

  void buildInfrastructure();
  void transferCoordinates();
  mfem::Array<int> mergedAttributes() const;

  const mfem::ParMesh& parent_mesh_;
  const mfem::ParGridFunction* coordinates_{};
  long parent_mesh_sequence_{};
  long parent_space_sequence_{};
  std::vector<int> mortar_attributes_;
  std::vector<int> nonmortar_attributes_;
  bool self_contact_{};
  int adjacency_exclusion_depth_{ 1 };
  MfemBridgeParameters parameters_;
  bool use_lor_{};
  int lor_factor_{};
  int maximum_order_{ native_high_order_max_order };
  int maximum_surface_order_{};

  std::unique_ptr<mfem::ParSubMesh> submesh_;
  std::unique_ptr<mfem::FiniteElementCollection> submesh_coordinate_collection_;
  std::unique_ptr<mfem::ParFiniteElementSpace> submesh_coordinate_space_;
  std::unique_ptr<mfem::ParGridFunction> submesh_coordinates_;
  std::unique_ptr<mfem::ParMesh> lor_mesh_;
  std::unique_ptr<mfem::InterpolationGridTransfer> high_order_to_lor_transfer_;
  mfem::ParFiniteElementSpace* surface_coordinate_space_{};
  mfem::ParGridFunction* surface_coordinates_{};
  std::unique_ptr<redecomp::RedecompMesh> redecomp_mesh_;
  std::unique_ptr<mfem::FiniteElementSpace> redecomp_coordinate_space_;
  std::unique_ptr<mfem::GridFunction> redecomp_coordinates_;
  std::unique_ptr<redecomp::RedecompTransfer> redecomp_transfer_;
  mfem::Array<int> submesh_to_parent_vdofs_;
  std::unique_ptr<ContactRestrictionOperator> coordinate_restriction_;
  std::unique_ptr<mfem::FiniteElementCollection> submesh_scalar_collection_;
  std::unique_ptr<mfem::ParFiniteElementSpace> submesh_scalar_space_;
  std::unique_ptr<mfem::FiniteElementCollection> lor_scalar_collection_;
  std::unique_ptr<mfem::ParFiniteElementSpace> surface_scalar_space_owner_;
  mfem::ParFiniteElementSpace* surface_scalar_space_{};
  std::unique_ptr<mfem::InterpolationGridTransfer> scalar_high_order_to_lor_transfer_;
  std::unique_ptr<mfem::FiniteElementSpace> redecomp_scalar_space_;
  std::unique_ptr<redecomp::RedecompTransfer> scalar_redecomp_transfer_;
  std::unique_ptr<SurfaceRestrictionOperator> scalar_restriction_;
  ContactDomain domain_;

  GeometryVersion geometry_version_;
  InteractionVersion interaction_version_;
  bool interactions_valid_{};
};

}  // namespace tribol::future

#endif
