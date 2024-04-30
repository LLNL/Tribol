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

namespace tribol
{

enum class MemorySpace
{
  Dynamic,
#ifdef TRIBOL_USE_UMPIRE
  Host,
  Device,
  Unified
#endif
};

enum class ExecutionMode
{
   Sequential,
#ifdef TRIBOL_USE_RAJA
#ifdef TRIBOL_USE_OPENMP
   SharedMemParallel,
#endif
#ifdef TRIBOL_USE_CUDA
   Cuda,
#endif
#ifdef TRIBOL_USE_HIP
   Hip,
#endif
   Dynamic
#endif
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
#ifdef TRIBOL_USE_UMPIRE
  return axom::getUmpireResourceAllocatorID(toUmpireMemoryType(mem_space));
#else
  return 0;
#endif
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

} // namespace tribol


#endif /* SRC_COMMON_EXECMODEL_HPP_ */
