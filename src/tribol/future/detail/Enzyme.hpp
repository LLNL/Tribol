#ifndef TRIBOL_FUTURE_DETAIL_ENZYME_HPP_
#define TRIBOL_FUTURE_DETAIL_ENZYME_HPP_

#include "tribol/future/Config.hpp"

#ifdef TRIBOL_USE_ENZYME
#ifdef MFEM_USE_ENZYME
#include "mfem/general/enzyme.hpp"
#else
extern int enzyme_dup;
extern int enzyme_dupnoneed;
extern int enzyme_out;
extern int enzyme_const;

#if defined( MFEM_USE_CUDA ) || defined( MFEM_USE_HIP )
extern __device__ int enzyme_dup;
extern __device__ int enzyme_dupnoneed;
extern __device__ int enzyme_out;
extern __device__ int enzyme_const;
#endif

template <typename Return, typename... Arguments>
MFEM_HOST_DEVICE Return __enzyme_autodiff( Arguments... );

template <typename Return, typename... Arguments>
MFEM_HOST_DEVICE Return __enzyme_fwddiff( Arguments... );
#endif

#if !defined( TRIBOL_USE_HOST ) && !defined( TRIBOL_DEVICE_CODE )
extern "C" {
extern int tribol_future_host_enzyme_const asm( "enzyme_const" );
extern int tribol_future_host_enzyme_dup asm( "enzyme_dup" );
extern int tribol_future_host_enzyme_out asm( "enzyme_out" );
}
#define TRIBOL_FUTURE_ENZYME_CONST tribol_future_host_enzyme_const
#define TRIBOL_FUTURE_ENZYME_DUP tribol_future_host_enzyme_dup
#define TRIBOL_FUTURE_ENZYME_OUT tribol_future_host_enzyme_out
#else
#define TRIBOL_FUTURE_ENZYME_CONST enzyme_const
#define TRIBOL_FUTURE_ENZYME_DUP enzyme_dup
#define TRIBOL_FUTURE_ENZYME_OUT enzyme_out
#endif
#endif

#endif
