#include "tribol/future/TribolFuture.hpp"

#include <mpi.h>

#include <algorithm>

namespace {

using namespace tribol::future;

using BuiltInPenalty = EnergyMortar<NodalConstraints, QuadraticPenaltyLaw, EnzymeLinearization>;
using BuiltInContact = MfemContactOperator<BuiltInPenalty, LowOrderRefinedSurface, SequentialExecution, BvhSearch>;

using HostPressure = EnergyMortar<NodalConstraints, ExternalPressureLaw, EnzymeLinearization>;
using HostPressureContact = MfemContactOperator<HostPressure, LowOrderRefinedSurface, SequentialExecution, BvhSearch>;

using Multiplier = EnergyMortar<NodalConstraints, NodalLagrangeMultiplier, EnzymeLinearization>;
using MultiplierContact = MfemContactOperator<Multiplier, LowOrderRefinedSurface, SequentialExecution, BvhSearch>;

[[maybe_unused]] void addBuiltInPenaltyResidual( mfem::ParMesh& mesh, mfem::ParGridFunction& coordinates,
                                                 mfem::Vector& residual )
{
  BuiltInContact::Parameters parameters;
  parameters.algorithm.enforcement.stiffness = 1.0e4;
  parameters.algorithm.quadrature_points = 3;
  BuiltInContact contact( mesh, coordinates, PairedSurfaces{ { 1 }, { 2 } }, parameters );
  contact.updateInteractions();
  contact.addResidual( residual );

  contact.updateGeometry( coordinates );
  contact.addResidual( residual );

  auto tangent = contact.linearize();
  mfem::Vector direction( residual.Size() );
  mfem::Vector tangent_action( residual.Size() );
  direction = 0.0;
  tangent.dforce_dx->Mult( direction, tangent_action );
}

[[maybe_unused]] void addHostPressureResidual( mfem::ParMesh& mesh, mfem::ParGridFunction& coordinates,
                                               mfem::Vector& residual )
{
  HostPressureContact contact( mesh, coordinates, PairedSurfaces{ { 1 }, { 2 } } );
  contact.updateInteractions();
  const auto nodal = contact.evaluateNodalKinematics();

  mfem::Vector potential( nodal.gap.Size() );
  mfem::Vector pressure( nodal.gap.Size() );
  mfem::Vector tangent( nodal.gap.Size() );
  const Real* gap = nodal.gap.HostRead();
  Real* psi = potential.HostWrite();
  Real* traction = pressure.HostWrite();
  Real* stiffness = tangent.HostWrite();
  for ( int dof = 0; dof < nodal.gap.Size(); ++dof ) {
    const Real active_gap = std::min( gap[dof], 0.0 );
    psi[dof] = 0.5e4 * active_gap * active_gap;
    traction[dof] = 1.0e4 * active_gap;
    stiffness[dof] = gap[dof] < 0.0 ? 1.0e4 : 0.0;
  }
  ExternalPressureData pressure_data{ potential, pressure, tangent };
  contact.addResidual( pressure, residual );
  auto linearization = contact.linearize( pressure_data );
  (void)linearization;
}

[[maybe_unused]] void addMultiplierResidual( mfem::ParMesh& mesh, mfem::ParGridFunction& coordinates,
                                             mfem::Vector& residual )
{
  MultiplierContact contact( mesh, coordinates, PairedSurfaces{ { 1 }, { 2 } } );
  contact.updateInteractions();

  mfem::Vector multiplier( contact.bridge().scalarRestriction().Width() );
  multiplier = 0.0;
  const auto result = contact.addResidual( multiplier, residual );
  const auto linearization = contact.linearize( multiplier );

  mfem::Vector coordinate_direction( residual.Size() );
  mfem::Vector force_direction( residual.Size() );
  mfem::Vector multiplier_direction( multiplier.Size() );
  mfem::Vector constraint_direction( multiplier.Size() );
  coordinate_direction = 0.0;
  multiplier_direction = 0.0;
  linearization.dforce_dx->Mult( coordinate_direction, force_direction );
  linearization.dforce_dlambda->Mult( multiplier_direction, force_direction );
  linearization.dgap_dx->Mult( coordinate_direction, constraint_direction );
  (void)result;
}

}  // namespace

int main( int argc, char* argv[] )
{
  MPI_Init( &argc, &argv );
  MPI_Finalize();
  return 0;
}
