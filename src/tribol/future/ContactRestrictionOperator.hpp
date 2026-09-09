#ifndef TRIBOL_FUTURE_CONTACTRESTRICTIONOPERATOR_HPP_
#define TRIBOL_FUTURE_CONTACTRESTRICTIONOPERATOR_HPP_

#include "tribol/future/Config.hpp"

#include "redecomp/RedecompTransfer.hpp"

namespace tribol::future {

class ContactRestrictionOperator final : public mfem::Operator {
 public:
  ContactRestrictionOperator( const mfem::ParFiniteElementSpace& parent_space,
                              mfem::ParFiniteElementSpace& submesh_space, mfem::ParFiniteElementSpace& surface_space,
                              mfem::FiniteElementSpace& redecomp_space, const mfem::Array<int>& submesh_to_parent_vdofs,
                              const mfem::Operator* high_order_to_lor,
                              const redecomp::RedecompTransfer& redecomp_transfer, bool use_device );

  void Mult( const mfem::Vector& parent_true_dofs, mfem::Vector& contact_dofs ) const override;
  void MultTranspose( const mfem::Vector& contact_dual, mfem::Vector& parent_true_dual ) const override;

 private:
  void zeroNonOwnedSurfaceDofs() const;

  const mfem::ParFiniteElementSpace& parent_space_;
  mfem::ParFiniteElementSpace& submesh_space_;
  mfem::ParFiniteElementSpace& surface_space_;
  mfem::FiniteElementSpace& redecomp_space_;
  const mfem::Array<int>& submesh_to_parent_vdofs_;
  const mfem::Operator* high_order_to_lor_{};
  const redecomp::RedecompTransfer& redecomp_transfer_;
  bool use_device_{};

  mutable mfem::Vector parent_local_;
  mutable mfem::Vector submesh_values_;
  mutable mfem::Vector surface_values_;
  mutable mfem::Vector redecomp_values_;
};

}  // namespace tribol::future

#endif
