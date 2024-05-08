// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Tribol Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (MIT)
 
#ifndef SRC_COMMON_LOOPEXEC_HPP_
#define SRC_COMMON_LOOPEXEC_HPP_

// C++ includes
#include <type_traits>

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
  struct forAllType {};

  template <ExecutionMode EXEC, typename BODY, bool ASYNC, int BLOCK_SIZE>
  void forAllImpl(forAllType<EXEC>, IndexT, BODY&&)
  {
    SLIC_ERROR_ROOT("forAllExec not defined for the given ExecutionMode.");
  }

#ifndef TRIBOL_USE_RAJA
// ExecutionMode::Dynamic maps to ExecutionMode::Sequential without RAJA
  template <typename BODY, bool ASYNC, int BLOCK_SIZE>
  void forAllImpl(forAllType<ExecutionMode::Dynamic>, IndexT N, BODY&& body)
  {
    forAllImpl(forAllType<ExecutionMode::Sequential>(), N, std::move(body));
  }
#endif

  template <typename BODY, bool ASYNC, int BLOCK_SIZE>
  void forAllImpl(forAllType<ExecutionMode::Sequential>, IndexT N, BODY&& body)
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
  template <typename BODY, bool ASYNC, int BLOCK_SIZE>
  typename std::enable_if_t<ASYNC> forAllCudaImpl(IndexT N, BODY&& body)
  {
    RAJA::forall<RAJA::cuda_exec_async<BLOCK_SIZE>>(RAJA::RangeSegment(0, N), std::move(body));
  }

  template <typename BODY, bool ASYNC, int BLOCK_SIZE>
  typename std::enable_if_t<!ASYNC> forAllCudaImpl(IndexT N, BODY&& body)
  {
    RAJA::forall<RAJA::cuda_exec<BLOCK_SIZE>>(RAJA::RangeSegment(0, N), std::move(body));
  }

  template <typename BODY, bool ASYNC, int BLOCK_SIZE>
  void forAllImpl(forAllType<ExecutionMode::Cuda>, IndexT N, BODY&& body)
  {
    forAllCudaImpl<BODY, ASYNC, BLOCK_SIZE>(N, std::move(body));
  }
#endif

#ifdef TRIBOL_USE_HIP
  template <typename BODY, bool ASYNC, int BLOCK_SIZE>
  typename std::enable_if_t<ASYNC> forAllHipImpl(IndexT N, BODY&& body)
  {
    RAJA::forall<RAJA::hip_exec_async<BLOCK_SIZE>>(RAJA::RangeSegment(0, N), std::move(body));
  }

  template <typename BODY, bool ASYNC, int BLOCK_SIZE>
  typename std::enable_if_t<!ASYNC> forAllHipImpl(IndexT N, BODY&& body)
  {
    RAJA::forall<RAJA::hip_exec<BLOCK_SIZE>>(RAJA::RangeSegment(0, N), std::move(body));
  }

  template <typename BODY, bool ASYNC, int BLOCK_SIZE>
  void forAllImpl(forAllType<ExecutionMode::Hip>, IndexT N, BODY&& body)
  {
    forAllHipImpl<BODY, ASYNC, BLOCK_SIZE>(N, std::move(body));
  }
#endif

#ifdef TRIBOL_USE_OPENMP
  template <typename BODY, bool ASYNC, int BLOCK_SIZE>
  void forAllImpl(forAllType<ExecutionMode::OpenMP>, IndexT N, BODY&& body)
  {
    RAJA::forall<RAJA::omp_parallel_for_exec>(RAJA::RangeSegment(0, N), std::move(body));
  }
#endif
}

#define TRIBOL_BLOCK_SIZE 256
template <ExecutionMode EXEC, typename BODY, bool ASYNC = false, int BLOCK_SIZE = TRIBOL_BLOCK_SIZE>
void forAllExec(IndexT N, BODY&& body)
{
  detail::forAllImpl<BODY, ASYNC, BLOCK_SIZE>(detail::forAllType<EXEC>(), N, std::move(body));
}

template <typename BODY, bool ASYNC = false, int BLOCK_SIZE = TRIBOL_BLOCK_SIZE>
void forAllExec(ExecutionMode exec_mode, IndexT N, BODY&& body)
{
  switch (exec_mode)
  {
    case ExecutionMode::Sequential:
      return detail::forAllImpl<BODY, ASYNC, BLOCK_SIZE>(
        detail::forAllType<ExecutionMode::Sequential>(), N, std::move(body));
#ifdef TRIBOL_USE_RAJA
#ifdef TRIBOL_USE_OPENMP
    case ExecutionMode::OpenMP:
      return detail::forAllImpl<BODY, ASYNC, BLOCK_SIZE>(
        detail::forAllType<ExecutionMode::OpenMP>(), N, std::move(body));
#endif
#ifdef TRIBOL_USE_CUDA
    case ExecutionMode::Cuda:
      return detail::forAllImpl<BODY, ASYNC, BLOCK_SIZE>(
        detail::forAllType<ExecutionMode::Cuda>(), N, std::move(body));
#endif
#ifdef TRIBOL_USE_HIP
    case ExecutionMode::Hip:
      return detail::forAllImpl<BODY, ASYNC, BLOCK_SIZE>(
        detail::forAllType<ExecutionMode::Hip>(), N, std::move(body));
#endif
#endif
    default:
      SLIC_ERROR_ROOT("Unsupported execution mode in a forAllExec loop.");
      return;
  }
}

} // namespace tribol


#endif /* SRC_COMMON_LOOPEXEC_HPP_ */
