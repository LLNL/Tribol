// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Tribol Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (MIT)
 
#ifndef SRC_COMMON_LOOPEXEC_HPP_
#define SRC_COMMON_LOOPEXEC_HPP_

// Tribol includes
#include "tribol/common/ExecModel.hpp"

// Axom includes
#include "axom/slic.hpp"

// RAJA includes
#ifdef TRIBOL_USE_RAJA
#include "RAJA/RAJA.hpp"
#endif

namespace tribol
{

// Check Tribol has RAJA if we are using CUDA, HIP, or OpenMP.
#if defined(TRIBOL_USE_CUDA) && !defined(TRIBOL_USE_RAJA)
#error "RAJA is required for CUDA support in tribol."
#endif 

#if defined(TRIBOL_USE_HIP) && !defined(TRIBOL_USE_RAJA)
#error "RAJA is required for HIP support in tribol."
#endif 

#if defined(TRIBOL_USE_OPENMP) && !defined(TRIBOL_USE_RAJA)
#error "RAJA is required for OpenMP support in tribol."
#endif 

// Check for compatibility with RAJA build configuration.
// RAJA_ENABLE_CUDA, RAJA_ENABLE_HIP, and RAJA_ENABLE_OPENMP are defined in RAJA/config.hpp.
#if defined(TRIBOL_USE_CUDA) && !defined(RAJA_ENABLE_CUDA)
#error "ENABLE_CUDA was specified for tribol, but RAJA was built without RAJA_ENABLE_CUDA"
#endif 

#if defined(TRIBOL_USE_HIP) && !defined(RAJA_ENABLE_HIP)
#error "ENABLE_HIP was specified for tribol, but RAJA was built without RAJA_ENABLE_HIP"
#endif 

#if defined(TRIBOL_USE_OPENMP) && !defined(RAJA_ENABLE_OPENMP)
#error "ENABLE_OPENMP was specified for tribol, but RAJA was built without RAJA_ENABLE_OPENMP"
#endif 

namespace detail
{
  // SFINAE type for choosing correct RAJA::forall policy
  template <ExecutionMode T>
  struct forallType {};

#ifdef TRIBOL_USE_RAJA
  template <ExecutionMode EXEC, typename BODY>
  void forAllImpl(forallType<EXEC>, const IndexT, BODY&&)
  {
    SLIC_ERROR_ROOT("forAllExec not defined for the given ExecutionMode.");
  }
#endif

#ifndef TRIBOL_USE_RAJA
// ExecutionMode::Dynamic maps to ExecutionMode::Sequential without RAJA
  template <typename BODY>
  void forAllImpl(forallType<ExecutionMode::Dynamic>, const IndexT N, BODY&& body)
  {
    forAllImpl(forallType<ExecutionMode::Sequential>(), N, std::move(body));
  }
#endif

  template <typename BODY>
  void forAllImpl(forallType<ExecutionMode::Sequential>, const IndexT N, BODY&& body)
  {
#ifdef TRIBOL_USE_RAJA
    RAJA::forall<RAJA::loop_exec>(RAJA::RangeSegment(0, N), std::move(body));
#else
    for (IndexT i{0}; i < N; ++i)
    {
      body(i);
    }
#endif
  }

#ifdef TRIBOL_USE_CUDA
#define TRIBOL_CUDA_BLOCK_SIZE 256
  template <const int BLOCK_SIZE=TRIBOL_CUDA_BLOCK_SIZE, typename BODY>
  void forAllCudaImpl(const IndexT N, BODY&& body)
  {
    RAJA::forall<RAJA::cuda_exec<BLOCK_SIZE>>(RAJA::RangeSegment(0, N), std::move(body));
  }

  template <typename BODY>
  void forAllImpl(forallType<ExecutionMode::Cuda>, const IndexT N, BODY&& body)
  {
    forAllCudaImpl<TRIBOL_CUDA_BLOCK_SIZE>(N, std::move(body));
  }
#endif

#ifdef TRIBOL_USE_HIP
#define TRIBOL_HIP_BLOCK_SIZE 256
  template <typename BODY>
  void forAllImpl(forallType<ExecutionMode::Hip>, const IndexT N, BODY&& body)
  {
    forAllHipImpl<TRIBOL_HIP_BLOCK_SIZE>(N, std::move(body));
  }

  template <const int BLOCK_SIZE=TRIBOL_HIP_BLOCK_SIZE, typename BODY>
  void forAllHipImpl(const IndexT N, BODY&& body)
  {
    RAJA::forall<RAJA::hip_exec<BLOCK_SIZE>>(RAJA::RangeSegment(0, N), std::move(body));
  }
#endif

#ifdef TRIBOL_USE_OPENMP
  template <typename BODY>
  void forAllImpl(forallType<ExecutionMode::SharedMemParallel>, const IndexT N, BODY&& body)
  {
    RAJA::forall<RAJA::omp_parallel_for_exec>(RAJA::RangeSegment(0, N), std::move(body));
  }
#endif
}

template <ExecutionMode EXEC, typename BODY>
void forAllExec(const IndexT N, BODY&& body)
{
  detail::forAllImpl(detail::forallType<EXEC>(), N, std::move(body));
}

} // namespace tribol


#endif /* SRC_COMMON_LOOPEXEC_HPP_ */
