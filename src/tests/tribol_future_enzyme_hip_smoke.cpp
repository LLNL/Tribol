#include "gtest/gtest.h"

#include "tribol/future/detail/Enzyme.hpp"

namespace tribol::future {
namespace {

MFEM_HOST_DEVICE void localEnergy( const Real* coordinates, Real* energy )
{
  *energy =
      0.5 * coordinates[0] * coordinates[0] + coordinates[0] * coordinates[1] + 1.5 * coordinates[1] * coordinates[1];
}

MFEM_HOST_DEVICE void localGradient( const Real* coordinates, Real* gradient )
{
  Real coordinate_shadow[2] = { 0.0, 0.0 };
  Real energy{};
  Real energy_shadow{ 1.0 };
  __enzyme_autodiff<void>( (void*)localEnergy, TRIBOL_FUTURE_ENZYME_DUP, coordinates, coordinate_shadow,
                           TRIBOL_FUTURE_ENZYME_DUP, &energy, &energy_shadow );
  gradient[0] = coordinate_shadow[0];
  gradient[1] = coordinate_shadow[1];
}

TEST( FutureEnzymeHipSmoke, DeviceGradientAndHessianVectorProduct )
{
  mfem::Vector coordinates( 2 );
  mfem::Vector direction( 2 );
  mfem::Vector gradient( 2 );
  mfem::Vector hessian_vector( 2 );
  coordinates.UseDevice( true );
  direction.UseDevice( true );
  gradient.UseDevice( true );
  hessian_vector.UseDevice( true );
  coordinates[0] = 2.0;
  coordinates[1] = -1.0;
  direction[0] = 0.25;
  direction[1] = -0.5;

  const Real* x = coordinates.Read( true );
  const Real* dx = direction.Read( true );
  Real* grad = gradient.Write( true );
  Real* hess_vec = hessian_vector.Write( true );
  mfem::forall( 1, [=] MFEM_HOST_DEVICE( int ) {
    localGradient( x, grad );
    __enzyme_fwddiff<void>( (void*)localGradient, TRIBOL_FUTURE_ENZYME_DUP, x, dx, TRIBOL_FUTURE_ENZYME_DUP, grad,
                            hess_vec );
  } );
  MFEM_DEVICE_SYNC;

  EXPECT_NEAR( gradient.HostRead()[0], 1.0, 1.0e-12 );
  EXPECT_NEAR( gradient.HostRead()[1], -1.0, 1.0e-12 );
  EXPECT_NEAR( hessian_vector.HostRead()[0], -0.25, 1.0e-12 );
  EXPECT_NEAR( hessian_vector.HostRead()[1], -1.25, 1.0e-12 );
}

}  // namespace
}  // namespace tribol::future
