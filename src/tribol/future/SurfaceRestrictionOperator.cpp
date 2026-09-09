#include "tribol/future/SurfaceRestrictionOperator.hpp"

#include <stdexcept>

namespace tribol::future {

SurfaceRestrictionOperator::SurfaceRestrictionOperator( mfem::ParFiniteElementSpace& coarse_space,
                                                        mfem::ParFiniteElementSpace& surface_space,
                                                        mfem::FiniteElementSpace& redecomp_space,
                                                        const mfem::Operator* coarse_to_surface,
                                                        const redecomp::RedecompTransfer& redecomp_transfer,
                                                        bool use_device )
    : mfem::Operator( redecomp_space.GetVSize(), coarse_space.GetTrueVSize() ),
      coarse_space_( coarse_space ),
      surface_space_( surface_space ),
      redecomp_space_( redecomp_space ),
      coarse_to_surface_( coarse_to_surface ),
      redecomp_transfer_( redecomp_transfer ),
      use_device_( use_device ),
      coarse_local_( coarse_space.GetVSize() ),
      surface_local_( surface_space.GetVSize() ),
      redecomp_local_( redecomp_space.GetVSize() )
{
  coarse_local_.UseDevice( use_device_ );
  surface_local_.UseDevice( use_device_ );
  redecomp_local_.UseDevice( use_device_ );
}

void SurfaceRestrictionOperator::Mult( const mfem::Vector& surface_true_dofs, mfem::Vector& contact_dofs ) const
{
  if ( surface_true_dofs.Size() != Width() || contact_dofs.Size() != Height() ) {
    throw std::invalid_argument( "SurfaceRestrictionOperator::Mult received an incorrectly sized vector." );
  }
  coarse_space_.GetProlongationMatrix()->Mult( surface_true_dofs, coarse_local_ );
  if ( coarse_to_surface_ != nullptr ) {
    coarse_to_surface_->Mult( coarse_local_, surface_local_ );
  } else {
    surface_local_ = coarse_local_;
  }
  mfem::ParGridFunction surface_field( &surface_space_, surface_local_ );
  mfem::GridFunction contact_field( &redecomp_space_, contact_dofs );
  redecomp_transfer_.TransferToSerial( surface_field, contact_field );
  contact_dofs.SyncMemory( contact_field );
}

void SurfaceRestrictionOperator::MultTranspose( const mfem::Vector& contact_dual,
                                                mfem::Vector& surface_true_dual ) const
{
  if ( contact_dual.Size() != Height() || surface_true_dual.Size() != Width() ) {
    throw std::invalid_argument( "SurfaceRestrictionOperator::MultTranspose received an incorrectly sized vector." );
  }
  redecomp_local_ = contact_dual;
  surface_local_ = 0.0;
  mfem::GridFunction contact_field( &redecomp_space_, redecomp_local_ );
  mfem::ParGridFunction surface_field( &surface_space_, surface_local_ );
  redecomp_transfer_.TransferToParallel( contact_field, surface_field );
  surface_local_.SyncMemory( surface_field );
  zeroNonOwnedSurfaceDofs();
  if ( coarse_to_surface_ != nullptr ) {
    coarse_local_ = 0.0;
    coarse_to_surface_->MultTranspose( surface_local_, coarse_local_ );
  } else {
    coarse_local_ = surface_local_;
  }
  coarse_space_.GetProlongationMatrix()->MultTranspose( coarse_local_, surface_true_dual );
}

void SurfaceRestrictionOperator::zeroNonOwnedSurfaceDofs() const
{
  const auto* prolongation = surface_space_.Dof_TrueDof_Matrix();
  if ( prolongation == nullptr ) {
    return;
  }
  const int* diagonal_offsets =
      mfem::Read( prolongation->GetDiagMemoryI(), surface_space_.GetVSize() + 1, use_device_ );
  Real* values = surface_local_.ReadWrite( use_device_ );
  mfem::forall_switch( use_device_, surface_local_.Size(), [=] MFEM_HOST_DEVICE( int index ) {
    if ( diagonal_offsets[index + 1] == diagonal_offsets[index] ) {
      values[index] = 0.0;
    }
  } );
}

}  // namespace tribol::future
