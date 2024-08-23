// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Tribol Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (MIT)
 
#ifndef SRC_COMMON_EXECMODEL_HPP_
#define SRC_COMMON_EXECMODEL_HPP_

// Tribol includes
#include "tribol/common/BasicTypes.hpp"

// Axom includes
#include "axom/core/memory_management.hpp"
#include "axom/slic.hpp"

namespace tribol
{

enum class MemorySpace
{
  Dynamic,
  Host,
#ifdef TRIBOL_USE_UMPIRE
  Device,
  Unified
#endif
};

enum class ExecutionMode
{
  Sequential,
#ifdef TRIBOL_USE_RAJA
#ifdef TRIBOL_USE_OPENMP
  OpenMP,
#endif
#ifdef TRIBOL_USE_CUDA
  Cuda,
#endif
#ifdef TRIBOL_USE_HIP
  Hip,
#endif
#endif
  Dynamic
};

template <MemorySpace MSPACE>
struct toAxomMemorySpace
{
  static constexpr axom::MemorySpace value = axom::MemorySpace::Dynamic;
};

#ifdef TRIBOL_USE_UMPIRE

template <>
struct toAxomMemorySpace<MemorySpace::Host>
{
  static constexpr axom::MemorySpace value = axom::MemorySpace::Host;
};

template <>
struct toAxomMemorySpace<MemorySpace::Device>
{
  static constexpr axom::MemorySpace value = axom::MemorySpace::Device;
};

template <>
struct toAxomMemorySpace<MemorySpace::Unified>
{
  static constexpr axom::MemorySpace value = axom::MemorySpace::Unified;
};

#endif

// Waiting for C++ 17...
// template <MemorySpace MSPACE>
// inline constexpr axom::MemorySpace axomMemorySpaceV = toAxomMemorySpace<MSPACE>::value;

#ifdef TRIBOL_USE_UMPIRE

inline umpire::resource::MemoryResourceType toUmpireMemoryType(MemorySpace mem_space)
{
  switch (mem_space)
  {
    case MemorySpace::Host:
      return umpire::resource::MemoryResourceType::Host;
    case MemorySpace::Device:
      return umpire::resource::MemoryResourceType::Device;
    case MemorySpace::Unified:
      return umpire::resource::MemoryResourceType::Unified;
    default:
      return umpire::resource::MemoryResourceType::Unknown;
  }
};

#endif

inline int getResourceAllocatorID(MemorySpace mem_space)
{
  int allocator_id = 0;
#ifdef TRIBOL_USE_UMPIRE
  if (mem_space == MemorySpace::Dynamic)
  {
    allocator_id = axom::getDefaultAllocatorID();
  }
  else
  {
    allocator_id = axom::getUmpireResourceAllocatorID(toUmpireMemoryType(mem_space));
  }
#else
  TRIBOL_UNUSED_VAR(mem_space);
#endif
  return allocator_id;
}

template <MemorySpace MSPACE>
struct toExecutionMode
{
#ifdef TRIBOL_USE_RAJA
  static constexpr ExecutionMode value = ExecutionMode::Sequential;
#else
  static constexpr ExecutionMode value = ExecutionMode::Sequential;
#endif
};

#ifdef TRIBOL_USE_CUDA
template <>
struct toExecutionMode<MemorySpace::Device>
{
  static constexpr ExecutionMode value = ExecutionMode::Cuda;
};
#endif

#ifdef TRIBOL_USE_HIP
template <>
struct toExecutionMode<MemorySpace::Device>
{
  static constexpr ExecutionMode value = ExecutionMode::Hip;
};
#endif

// Waiting for C++ 17...
// template <MemorySpace MSPACE>
// inline constexpr ExecutionMode ExecutionModeV = toExecutionMode<MSPACE>::value;

inline ExecutionMode getExecutionMode(MemorySpace mem_space)
{
  switch (mem_space)
  {
#ifdef TRIBOL_USE_UMPIRE
    case MemorySpace::Device:
  #ifdef TRIBOL_USE_RAJA
    #if defined(TRIBOL_USE_CUDA)
      return ExecutionMode::Cuda;
    #elif defined(TRIBOL_USE_HIP)
      return ExecutionMode::Hip;
    #else
      SLIC_WARNING_ROOT("No device execution mode available. "
        "Trying to run sequentially on host...");
      return ExecutionMode::Sequential;
    #endif
  #endif
    case MemorySpace::Unified:
#endif
    case MemorySpace::Dynamic:
#ifdef TRIBOL_USE_UMPIRE
      SLIC_DEBUG_ROOT("Dynamic or unified memory does not correspond to an execution space. "
        "Trying to run sequentially on host...");
    case MemorySpace::Host:
#endif
    default:
      return ExecutionMode::Sequential;
  }
}

inline ExecutionMode getExecutionMode(MemorySpace mem_space, ExecutionMode suggested_exec_mode)
{
  auto exec_mode = ExecutionMode::Sequential;
#ifdef TRIBOL_USE_RAJA
  switch (mem_space)
  {
    case MemorySpace::Dynamic:
  #ifdef TRIBOL_USE_UMPIRE
      // trust the user here...
      exec_mode = suggested_exec_mode;
      if (exec_mode == ExecutionMode::Dynamic)
      {
        SLIC_WARNING_ROOT("Dynamic execution with dynamic memory space. "
          "Assuming sequential execution on host.");
        exec_mode = ExecutionMode::Sequential;
      }
  #else
      // if we have RAJA but no Umpire, execute serially on host
      exec_mode = ExecutionMode::Sequential;
  #endif
      break;
  #ifdef TRIBOL_USE_UMPIRE
    case MemorySpace::Unified:
      // this should be able to run anywhere. let the user decide.
      exec_mode = suggested_exec_mode;
      if (exec_mode == ExecutionMode::Dynamic)
      {
    #if defined(TRIBOL_USE_CUDA)
        SLIC_INFO_ROOT("Dynamic execution with unified memory space. "
          "Assuming CUDA parallel execution.");
        exec_mode = ExecutionMode::Cuda;
    #elif defined(TRIBOL_USE_HIP)
        SLIC_INFO_ROOT("Dynamic execution with unified memory space. "
          "Assuming HIP parallel execution.");
        exec_mode = ExecutionMode::Hip;
    #elif defined(TRIBOL_USE_OPENMP)
        SLIC_INFO_ROOT("Dynamic execution with unified memory space. "
          "Assuming OpenMP parallel execution.");
        exec_mode = ExecutionMode::OpenMP;
    #else
        SLIC_INFO_ROOT("Dynamic execution with unified memory space. "
          "Assuming sequential execution.");
        exec_mode = ExecutionMode::Sequential;
    #endif
      }
      break;
  #endif
    case MemorySpace::Host:
      switch (suggested_exec_mode)
      {
        case ExecutionMode::Sequential:
  #ifdef TRIBOL_USE_OPENMP
        case ExecutionMode::OpenMP:
  #endif
          exec_mode = suggested_exec_mode;
          break;
        case ExecutionMode::Dynamic:
  #ifdef TRIBOL_USE_OPENMP
          SLIC_INFO_ROOT("Dynamic execution with a host memory space. "
            "Assuming OpenMP parallel execution.");
          exec_mode = ExecutionMode::OpenMP;
  #else
          SLIC_INFO_ROOT("Dynamic execution with a host memory space. "
            "Assuming sequential execution.");
          exec_mode = ExecutionMode::Sequential;
  #endif
          break;
        default:
          SLIC_WARNING_ROOT("Unsupported execution mode for host memory. "
            "Switching to sequential execution.");
          exec_mode = ExecutionMode::Sequential;
          break;
      }
      break;
  #ifdef TRIBOL_USE_UMPIRE
    case MemorySpace::Device:
      switch (suggested_exec_mode)
      {
    #ifdef TRIBOL_USE_CUDA
        case ExecutionMode::Cuda:
    #endif
    #ifdef TRIBOL_USE_HIP
        case ExecutionMode::Hip:
    #endif
          exec_mode = suggested_exec_mode;
          break;
        case ExecutionMode::Dynamic:
    #if defined(TRIBOL_USE_CUDA)
          exec_mode = ExecutionMode::Cuda;
          break;
    #elif defined(TRIBOL_USE_HIP)
          exec_mode = ExecutionMode::Hip;
          break;
    #endif
        default:
          SLIC_WARNING_ROOT("Unknown execution mode for device memory. "
            "Trying host sequential execution.");
          exec_mode = ExecutionMode::Sequential;
          break;
      }
  #endif
  }
#endif
  return exec_mode;
}

inline bool isOnDevice(ExecutionMode exec)
{
  switch (exec)
  {
#ifdef TRIBOL_USE_RAJA
  #if defined(TRIBOL_USE_CUDA)
    case ExecutionMode::Cuda:
      return true;
  #elif defined(TRIBOL_USE_HIP)
    case ExecutionMode::Hip:
      return true;
  #endif
  #ifdef TRIBOL_USE_OPENMP
    case ExecutionMode::OpenMP:
      return false;
  #endif
    case ExecutionMode::Dynamic:
      SLIC_ERROR_ROOT("Dynamic execution mode does not define a memory space location.");
      return false;
#endif
    case ExecutionMode::Sequential:
      return false;
    default:
      SLIC_ERROR_ROOT("Unknown execution mode.");
      return false;
  }
}

} // namespace tribol


#endif /* SRC_COMMON_EXECMODEL_HPP_ */
