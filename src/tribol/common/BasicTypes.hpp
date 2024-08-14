// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Tribol Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (MIT)

#ifndef TRIBOL_COMMON_BASICTYPES_HPP_
#define TRIBOL_COMMON_BASICTYPES_HPP_

// Tribol includes
#include "tribol/config.hpp"

// Axom includes
#include "axom/core/Types.hpp"

// MPI includes
#ifdef TRIBOL_USE_MPI
#include <mpi.h>
#endif

namespace tribol
{

#ifdef TRIBOL_USE_MPI

using CommT = MPI_Comm;
#define TRIBOL_COMM_WORLD MPI_COMM_WORLD
#define TRIBOL_COMM_NULL MPI_COMM_NULL

#else

using CommT = int;
#define TRIBOL_COMM_WORLD 0
#define TRIBOL_COMM_NULL -1

#endif

// match index type used in axom (since data is held in axom data structures)
using IndexT = axom::IndexType;

#ifdef TRIBOL_USE_SINGLE_PRECISION

#error "Tribol does not support single precision."
using RealT = float;

#else

using RealT = double;

#endif

#define TRIBOL_UNUSED_VAR AXOM_UNUSED_VAR
#define TRIBOL_UNUSED_PARAM AXOM_UNUSED_PARAM

// Execution space specifiers
#if defined(TRIBOL_USE_CUDA) || defined(TRIBOL_USE_HIP)
  #ifndef __device__
    #error "TRIBOL_USE_CUDA or TRIBOL_USE_HIP but __device__ is undefined.  Check include files"
  #endif
  #define TRIBOL_DEVICE __device__
  #define TRIBOL_HOST_DEVICE __host__ __device__
#else
  #define TRIBOL_DEVICE
  #define TRIBOL_HOST_DEVICE
#endif

// Define variable when loops are computed on host
#if !(defined(TRIBOL_USE_CUDA) || defined(TRIBOL_USE_HIP))
  #define TRIBOL_USE_HOST
#endif

} // namespace tribol

#endif /* TRIBOL_COMMON_BASICTYPES_HPP_ */
