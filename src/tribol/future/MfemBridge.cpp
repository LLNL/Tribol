#include "tribol/future/MfemBridge.hpp"

#include <algorithm>
#include <stdexcept>

namespace tribol::future {

void MfemBridge::updateGeometry( const mfem::ParGridFunction& coordinates )
{
  if ( coordinates.ParFESpace() != coordinates_->ParFESpace() || parent_mesh_.GetSequence() != parent_mesh_sequence_ ||
       coordinates.ParFESpace()->GetSequence() != parent_space_sequence_ ) {
    throw std::invalid_argument(
        "MfemBridge::updateGeometry requires the original finite-element space; reconstruct after AMR or "
        "repartitioning." );
  }
  coordinates_ = &coordinates;
  transferCoordinates();
  domain_.updateCoordinates( *redecomp_coordinates_ );
  ++geometry_version_;
}

void MfemBridge::rebuildGeometry( const mfem::ParGridFunction& coordinates )
{
  if ( coordinates.ParFESpace()->GetParMesh() != &parent_mesh_ || parent_mesh_.GetSequence() != parent_mesh_sequence_ ||
       coordinates.ParFESpace()->GetSequence() != parent_space_sequence_ ) {
    throw std::invalid_argument(
        "MfemBridge requires fixed parent topology and finite-element spaces; reconstruct after AMR or "
        "repartitioning." );
  }
  coordinates_ = &coordinates;
  buildInfrastructure();
  domain_.rebuild( *redecomp_mesh_, *redecomp_coordinates_, mortar_attributes_, nonmortar_attributes_, self_contact_ );
  ++geometry_version_;
  interactions_valid_ = false;
}

void MfemBridge::buildInfrastructure()
{
  coordinate_restriction_.reset();
  scalar_restriction_.reset();
  scalar_redecomp_transfer_.reset();
  redecomp_scalar_space_.reset();
  scalar_high_order_to_lor_transfer_.reset();
  surface_scalar_space_owner_.reset();
  lor_scalar_collection_.reset();
  submesh_scalar_space_.reset();
  submesh_scalar_collection_.reset();
  redecomp_transfer_.reset();
  redecomp_coordinates_.reset();
  redecomp_coordinate_space_.reset();
  redecomp_mesh_.reset();
  high_order_to_lor_transfer_.reset();
  lor_mesh_.reset();
  submesh_coordinates_.reset();
  submesh_coordinate_space_.reset();
  submesh_coordinate_collection_.reset();
  submesh_.reset();

  auto attributes = mergedAttributes();
  submesh_ = std::make_unique<mfem::ParSubMesh>( mfem::ParSubMesh::CreateFromBoundary( parent_mesh_, attributes ) );
  submesh_->EnsureNodes();

  const auto& parent_space = *coordinates_->ParFESpace();
  submesh_coordinate_collection_.reset( parent_space.FEColl()->Clone( parent_space.GetMaxElementOrder() ) );
  submesh_coordinate_space_ = std::make_unique<mfem::ParFiniteElementSpace>(
      submesh_.get(), submesh_coordinate_collection_.get(), parent_space.GetVDim(), parent_space.GetOrdering() );
  submesh_coordinates_ = std::make_unique<mfem::ParGridFunction>( submesh_coordinate_space_.get() );
  submesh_coordinates_->UseDevice( parameters_.use_device );
  mfem::ParSubMesh::Transfer( *coordinates_, *submesh_coordinates_ );

  submesh_scalar_collection_.reset( parent_space.FEColl()->Clone( parent_space.GetMaxElementOrder() ) );
  submesh_scalar_space_ =
      std::make_unique<mfem::ParFiniteElementSpace>( submesh_.get(), submesh_scalar_collection_.get(), 1 );
  surface_scalar_space_ = submesh_scalar_space_.get();

  surface_coordinate_space_ = submesh_coordinate_space_.get();
  surface_coordinates_ = submesh_coordinates_.get();
  if ( use_lor_ && parent_space.GetMaxElementOrder() > 1 ) {
    const int factor = lor_factor_ > 0 ? lor_factor_ : parent_space.GetMaxElementOrder();
    lor_mesh_ = std::make_unique<mfem::ParMesh>(
        mfem::ParMesh::MakeRefined( *submesh_, factor, mfem::BasisType::ClosedUniform ) );
    lor_mesh_->EnsureNodes();
    surface_coordinates_ = dynamic_cast<mfem::ParGridFunction*>( lor_mesh_->GetNodes() );
    if ( surface_coordinates_ == nullptr ) {
      throw std::runtime_error( "The LOR mesh did not create a parallel nodal grid function." );
    }
    surface_coordinate_space_ = surface_coordinates_->ParFESpace();
    high_order_to_lor_transfer_ =
        std::make_unique<mfem::InterpolationGridTransfer>( *submesh_coordinate_space_, *surface_coordinate_space_ );
    lor_scalar_collection_ = std::make_unique<mfem::H1_FECollection>( 1, lor_mesh_->Dimension() );
    surface_scalar_space_owner_ =
        std::make_unique<mfem::ParFiniteElementSpace>( lor_mesh_.get(), lor_scalar_collection_.get(), 1 );
    surface_scalar_space_ = surface_scalar_space_owner_.get();
    scalar_high_order_to_lor_transfer_ =
        std::make_unique<mfem::InterpolationGridTransfer>( *submesh_scalar_space_, *surface_scalar_space_ );
  } else if ( !use_lor_ && parent_space.GetMaxElementOrder() > maximum_order_ ) {
    throw std::invalid_argument( "Native high-order contact order exceeds the configured maximum." );
  }

  maximum_surface_order_ = surface_coordinate_space_->GetMaxElementOrder();
  MPI_Allreduce( MPI_IN_PLACE, &maximum_surface_order_, 1, MPI_INT, MPI_MAX, parent_mesh_.GetComm() );

  if ( parameters_.ghost_distance >= 0.0 ) {
    redecomp_mesh_ =
        std::make_unique<redecomp::RedecompMesh>( *surface_coordinate_space_->GetParMesh(), parameters_.ghost_distance,
                                                  redecomp::RedecompMesh::RCB, parameters_.redecomp_ranks );
  } else {
    redecomp_mesh_ = std::make_unique<redecomp::RedecompMesh>(
        *surface_coordinate_space_->GetParMesh(), redecomp::RedecompMesh::RCB, parameters_.redecomp_ranks );
  }
  redecomp_mesh_->EnsureNodes();
  redecomp_coordinate_space_ = std::make_unique<mfem::FiniteElementSpace>(
      redecomp_mesh_.get(), surface_coordinate_space_->FEColl(), surface_coordinate_space_->GetVDim(),
      surface_coordinate_space_->GetOrdering() );
  redecomp_coordinates_ = std::make_unique<mfem::GridFunction>( redecomp_coordinate_space_.get() );
  redecomp_coordinates_->UseDevice( parameters_.use_device );
  redecomp_transfer_ =
      std::make_unique<redecomp::RedecompTransfer>( *surface_coordinate_space_, *redecomp_coordinate_space_ );
  redecomp_scalar_space_ =
      std::make_unique<mfem::FiniteElementSpace>( redecomp_mesh_.get(), surface_scalar_space_->FEColl(), 1 );
  scalar_redecomp_transfer_ =
      std::make_unique<redecomp::RedecompTransfer>( *surface_scalar_space_, *redecomp_scalar_space_ );

  mfem::SubMeshUtils::BuildVdofToVdofMap( *submesh_coordinate_space_, parent_space, submesh_->GetFrom(),
                                          submesh_->GetParentElementIDMap(), submesh_to_parent_vdofs_ );
  submesh_to_parent_vdofs_.GetMemory().UseDevice( parameters_.use_device );

  transferCoordinates();
  coordinate_restriction_ = std::make_unique<ContactRestrictionOperator>(
      parent_space, *submesh_coordinate_space_, *surface_coordinate_space_, *redecomp_coordinate_space_,
      submesh_to_parent_vdofs_, high_order_to_lor_transfer_ ? &high_order_to_lor_transfer_->ForwardOperator() : nullptr,
      *redecomp_transfer_, parameters_.use_device );
  scalar_restriction_ = std::make_unique<SurfaceRestrictionOperator>(
      *submesh_scalar_space_, *surface_scalar_space_, *redecomp_scalar_space_,
      scalar_high_order_to_lor_transfer_ ? &scalar_high_order_to_lor_transfer_->ForwardOperator() : nullptr,
      *scalar_redecomp_transfer_, parameters_.use_device );
}

void MfemBridge::transferCoordinates()
{
  *submesh_coordinates_ = 0.0;
  mfem::ParSubMesh::Transfer( *coordinates_, *submesh_coordinates_ );
  if ( high_order_to_lor_transfer_ != nullptr ) {
    *surface_coordinates_ = 0.0;
    high_order_to_lor_transfer_->ForwardOperator().Mult( *submesh_coordinates_, *surface_coordinates_ );
  }
  *redecomp_coordinates_ = 0.0;
  redecomp_transfer_->TransferToSerial( *surface_coordinates_, *redecomp_coordinates_ );
}

mfem::Array<int> MfemBridge::mergedAttributes() const
{
  std::vector<int> attributes = mortar_attributes_;
  attributes.insert( attributes.end(), nonmortar_attributes_.begin(), nonmortar_attributes_.end() );
  std::sort( attributes.begin(), attributes.end() );
  attributes.erase( std::unique( attributes.begin(), attributes.end() ), attributes.end() );
  mfem::Array<int> result( static_cast<int>( attributes.size() ) );
  for ( int i = 0; i < result.Size(); ++i ) {
    result[i] = attributes[static_cast<std::size_t>( i )];
  }
  return result;
}

}  // namespace tribol::future
