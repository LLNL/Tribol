#include "tribol/future/ContactRestrictionOperator.hpp"

#include <stdexcept>

namespace tribol::future {

ContactRestrictionOperator::ContactRestrictionOperator(
    const mfem::ParFiniteElementSpace& parent_space, mfem::ParFiniteElementSpace& submesh_space,
    mfem::ParFiniteElementSpace& surface_space, mfem::FiniteElementSpace& redecomp_space,
    const mfem::Array<int>& submesh_to_parent_vdofs, const mfem::Operator* high_order_to_lor,
    const redecomp::RedecompTransfer& redecomp_transfer, bool use_device )
    : mfem::Operator( redecomp_space.GetVSize(), parent_space.GetTrueVSize() ),
      parent_space_( parent_space ),
      submesh_space_( submesh_space ),
      surface_space_( surface_space ),
      redecomp_space_( redecomp_space ),
      submesh_to_parent_vdofs_( submesh_to_parent_vdofs ),
      high_order_to_lor_( high_order_to_lor ),
      redecomp_transfer_( redecomp_transfer ),
      use_device_( use_device ),
      parent_local_( parent_space.GetVSize() ),
      submesh_values_( submesh_space.GetVSize() ),
      surface_values_( surface_space.GetVSize() ),
      redecomp_values_( redecomp_space.GetVSize() )
{
  if ( submesh_to_parent_vdofs_.Size() != submesh_space_.GetVSize() ) {
    throw std::invalid_argument( "The submesh-to-parent map does not match the submesh finite-element space." );
  }
  parent_local_.UseDevice( use_device_ );
  submesh_values_.UseDevice( use_device_ );
  surface_values_.UseDevice( use_device_ );
  redecomp_values_.UseDevice( use_device_ );
}

void ContactRestrictionOperator::Mult( const mfem::Vector& parent_true_dofs, mfem::Vector& contact_dofs ) const
{
  if ( parent_true_dofs.Size() != Width() || contact_dofs.Size() != Height() ) {
    throw std::invalid_argument( "ContactRestrictionOperator::Mult received an incorrectly sized vector." );
  }

  parent_space_.GetProlongationMatrix()->Mult( parent_true_dofs, parent_local_ );
  const Real* parent = parent_local_.Read( use_device_ );
  const int* map = submesh_to_parent_vdofs_.Read( use_device_ );
  Real* submesh = submesh_values_.Write( use_device_ );
  mfem::forall_switch( use_device_, submesh_values_.Size(), [=] MFEM_HOST_DEVICE( int index ) {
    const int encoded = map[index];
    const int parent_dof = encoded >= 0 ? encoded : -1 - encoded;
    submesh[index] = ( encoded >= 0 ? 1.0 : -1.0 ) * parent[parent_dof];
  } );

  if ( high_order_to_lor_ != nullptr ) {
    high_order_to_lor_->Mult( submesh_values_, surface_values_ );
  } else {
    surface_values_ = submesh_values_;
  }

  mfem::ParGridFunction surface_field( &surface_space_, surface_values_ );
  mfem::GridFunction contact_field( &redecomp_space_, redecomp_values_ );
  redecomp_transfer_.TransferToSerial( surface_field, contact_field );
  redecomp_values_.SyncMemory( contact_field );

  const int number_of_dofs = redecomp_space_.GetNDofs();
  const int dimension = redecomp_space_.GetVDim();
  const auto ordering = redecomp_space_.GetOrdering();
  const Real* native_values = redecomp_values_.Read( use_device_ );
  Real* structure_of_arrays = contact_dofs.Write( use_device_ );
  mfem::forall_switch( use_device_, contact_dofs.Size(), [=] MFEM_HOST_DEVICE( int index ) {
    const int component = index / number_of_dofs;
    const int dof = index - component * number_of_dofs;
    const int native_index =
        ordering == mfem::Ordering::byNODES ? component * number_of_dofs + dof : dof * dimension + component;
    structure_of_arrays[index] = native_values[native_index];
  } );
}

void ContactRestrictionOperator::MultTranspose( const mfem::Vector& contact_dual, mfem::Vector& parent_true_dual ) const
{
  if ( contact_dual.Size() != Height() || parent_true_dual.Size() != Width() ) {
    throw std::invalid_argument( "ContactRestrictionOperator::MultTranspose received an incorrectly sized vector." );
  }

  const int number_of_dofs = redecomp_space_.GetNDofs();
  const int dimension = redecomp_space_.GetVDim();
  const auto ordering = redecomp_space_.GetOrdering();
  const Real* structure_of_arrays = contact_dual.Read( use_device_ );
  Real* native_values = redecomp_values_.Write( use_device_ );
  mfem::forall_switch( use_device_, contact_dual.Size(), [=] MFEM_HOST_DEVICE( int index ) {
    const int component = index / number_of_dofs;
    const int dof = index - component * number_of_dofs;
    const int native_index =
        ordering == mfem::Ordering::byNODES ? component * number_of_dofs + dof : dof * dimension + component;
    native_values[native_index] = structure_of_arrays[index];
  } );
  surface_values_ = 0.0;
  mfem::GridFunction contact_field( &redecomp_space_, redecomp_values_ );
  mfem::ParGridFunction surface_field( &surface_space_, surface_values_ );
  redecomp_transfer_.TransferToParallel( contact_field, surface_field );
  surface_values_.SyncMemory( surface_field );
  zeroNonOwnedSurfaceDofs();

  if ( high_order_to_lor_ != nullptr ) {
    submesh_values_ = 0.0;
    high_order_to_lor_->MultTranspose( surface_values_, submesh_values_ );
  } else {
    submesh_values_ = surface_values_;
  }

  parent_local_ = 0.0;
  const Real* submesh = submesh_values_.Read( use_device_ );
  const int* map = submesh_to_parent_vdofs_.Read( use_device_ );
  Real* parent = parent_local_.ReadWrite( use_device_ );
  mfem::forall_switch( use_device_, submesh_values_.Size(), [=] MFEM_HOST_DEVICE( int index ) {
    const int encoded = map[index];
    const int parent_dof = encoded >= 0 ? encoded : -1 - encoded;
    AtomicAdd( parent[parent_dof], ( encoded >= 0 ? 1.0 : -1.0 ) * submesh[index] );
  } );

  parent_space_.GetProlongationMatrix()->MultTranspose( parent_local_, parent_true_dual );
}

void ContactRestrictionOperator::zeroNonOwnedSurfaceDofs() const
{
  const auto* prolongation = surface_space_.Dof_TrueDof_Matrix();
  if ( prolongation == nullptr ) {
    return;
  }
  const int* diagonal_offsets =
      mfem::Read( prolongation->GetDiagMemoryI(), surface_space_.GetVSize() + 1, use_device_ );
  Real* values = surface_values_.ReadWrite( use_device_ );
  mfem::forall_switch( use_device_, surface_values_.Size(), [=] MFEM_HOST_DEVICE( int index ) {
    if ( diagonal_offsets[index + 1] == diagonal_offsets[index] ) {
      values[index] = 0.0;
    }
  } );
}

}  // namespace tribol::future
