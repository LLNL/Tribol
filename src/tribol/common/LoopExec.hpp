// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Tribol Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (MIT)
 
#ifndef SRC_COMMON_LOOPEXEC_HPP_
#define SRC_COMMON_LOOPEXEC_HPP_

// Tribol includes
#include "tribol/types.hpp"
#include "tribol/common/Parameters.hpp"

// RAJA includes
#ifdef TRIBOL_USE_RAJA
#include "RAJA/RAJA.hpp"
#endif

namespace tribol
{


// Check for compatibility with RAJA build configuration.
// RAJA_ENABLE_CUDA and RAJA_ENABLE_OPENMP are defined in RAJA/config.hpp.
#if defined(TRIBOL_ENABLE_CUDA) && !defined(RAJA_ENABLE_CUDA)
#error "ENABLE_CUDA was specified for tribol, but RAJA was built without RAJA_ENABLE_CUDA"
#endif 

#if defined(TRIBOL_ENABLE_OPENMP) && !defined(RAJA_ENABLE_OPENMP)
#error "ENABLE_OPENMP was specified for tribol, but RAJA was built without RAJA_ENABLE_OPENMP"
#endif 

// This template executes a lambda executed in a sequential RAJA loop on the host
template <typename HBODY>
void RajaSeqExec(const int N, HBODY &&h_body)
{
   RAJA::forall<RAJA::loop_exec>(RAJA::RangeSegment(0,N), h_body);
}

#ifdef TRIBOL_ENABLE_CUDA
#define TRIBOL_CUDA_BLOCK_SIZE 256
// This template executes a lambda executed in a parallel RAJA loop on the device
template <const int BLOCK_SIZE=TRIBOL_CUDA_BLOCK_SIZE, typename DBODY>
void RajaCudaExec(const int N, DBODY &&d_body)
{
   // Note: a false argument is supplied to the cuda_exec policy so that 
   //       kernels will run synchronously.  Ohtherwise, the forall statement
   //       would need to be followed by RAJA::cuda::synchronize() to obtain
   //       similar behavior as the OpenMP forall.
   RAJA::forall<RAJA::cuda_exec<BLOCK_SIZE,false>>(RAJA::RangeSegment(0,N),d_body);
}
#endif

#ifdef TRIBOL_ENABLE_OPENMP
// This template wraps a lambda executed in a parallel RAJA loop on the host using OpenMP
template <typename HBODY>
void RajaOmpExec(const int N, HBODY &&h_body)
{
   RAJA::forall<RAJA::omp_parallel_for_exec>(RAJA::RangeSegment(0,N), h_body);
}
#endif

// This template function allows the loop execution type to be
// selected at run time.  The DBODY and HBODY paramters are
// lambdas to be run on the device and host, respectively.
template <typename DBODY, typename HBODY>
inline void ForallExec(const int N, DBODY&& d_body, HBODY&& h_body)
{
#ifdef TRIBOL_ENABLE_CUDA
   if (parameters_t::getInstance().exec_mode == CUDA_PARALLEL) { return RajaCudaExec(N, d_body); }
#endif
#ifdef TRIBOL_ENABLE_OPENMP
   if (parameters_t::getInstance().exec_mode == OPENMP_PARALLEL) { return RajaOmpExec(N, h_body); }
#endif
   return RajaSeqExec(N, h_body); // Default
}

/*!
   * \brief Perform a single loop on GPU or CPU with policy set by parameters.exec_mode 
   *
   * \param [in] i loop index variable
   * \param [in] N number of loop iterations, 0 <= i< N
   * \param [in] d_body lambda body to be executed on the device when parameters.exec_mode is CUDA_PARALLEL
   * \param [in] h_body lambda body to be executed on the host when parameters.exec_mode is OPENMP_PARALLEL or SEQUENTIAL
   *
   */
#define TRIBOL_FORALL(i,N,...)                         \
   ForallExec(N,                                       \
             [=] TRIBOL_DEVICE (int i) {__VA_ARGS__},  \
             [&] (int i) {__VA_ARGS__})

} // namespace tribol


#endif /* SRC_COMMON_LOOPEXEC_HPP_ */
