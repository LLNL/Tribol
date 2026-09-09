#ifndef TRIBOL_FUTURE_SURFACERESTRICTIONOPERATOR_HPP_
#define TRIBOL_FUTURE_SURFACERESTRICTIONOPERATOR_HPP_

#include "tribol/future/Config.hpp"

#include "redecomp/RedecompTransfer.hpp"

namespace tribol::future {

class SurfaceRestrictionOperator final : public mfem::Operator {
 public:
  SurfaceRestrictionOperator( mfem::ParFiniteElementSpace& coarse_space, mfem::ParFiniteElementSpace& surface_space,
                              mfem::FiniteElementSpace& redecomp_space, const mfem::Operator* coarse_to_surface,
                              const redecomp::RedecompTransfer& redecomp_transfer, bool use_device );

  void Mult( const mfem::Vector& surface_true_dofs, mfem::Vector& contact_dofs ) const override;
  void MultTranspose( const mfem::Vector& contact_dual, mfem::Vector& surface_true_dual ) const override;

 private:
  void zeroNonOwnedSurfaceDofs() const;

  mfem::ParFiniteElementSpace& coarse_space_;
  mfem::ParFiniteElementSpace& surface_space_;
  mfem::FiniteElementSpace& redecomp_space_;
  const mfem::Operator* coarse_to_surface_{};
  const redecomp::RedecompTransfer& redecomp_transfer_;
  bool use_device_{};
  mutable mfem::Vector coarse_local_;
  mutable mfem::Vector surface_local_;
  mutable mfem::Vector redecomp_local_;
};

}  // namespace tribol::future

#endif
