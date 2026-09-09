#ifndef TRIBOL_FUTURE_SURFACEDISCRETIZATIONS_HPP_
#define TRIBOL_FUTURE_SURFACEDISCRETIZATIONS_HPP_

#include "tribol/future/Config.hpp"

#include <stdexcept>

namespace tribol::future {

class MfemBridge;

struct LowOrderRefinedSurface {
  struct Parameters {
    int refinement_factor{ 0 };
  };

  static constexpr bool is_low_order_refined = true;

  static void validate( const mfem::ParGridFunction&, const Parameters& parameters )
  {
    if ( parameters.refinement_factor < 0 ) {
      throw std::invalid_argument( "The LOR refinement factor cannot be negative." );
    }
  }
};

template <int MaxOrder = native_high_order_max_order>
struct NativeHighOrderSurface {
  static_assert( MaxOrder > 0 );

  struct Parameters {};

  static constexpr bool is_low_order_refined = false;
  static constexpr int max_order = MaxOrder;

  static void validate( const mfem::ParGridFunction& coordinates, const Parameters& )
  {
    int order = coordinates.ParFESpace()->GetMaxElementOrder();
    MPI_Allreduce( MPI_IN_PLACE, &order, 1, MPI_INT, MPI_MAX, coordinates.ParFESpace()->GetComm() );
    if ( order > MaxOrder ) {
      throw std::invalid_argument( "The coordinate finite-element order exceeds NativeHighOrderSurface::MaxOrder." );
    }
  }
};

}  // namespace tribol::future

#endif
