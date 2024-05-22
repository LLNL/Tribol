// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Tribol Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (MIT)

#include <gtest/gtest.h>

#include "tribol/interface/tribol.hpp"


class ExecutionModeTest : public testing::TestWithParam<std::tuple<tribol::MemorySpace,    // given memory space
                                                                   tribol::ExecutionMode,  // given execution mode
                                                                   tribol::ExecutionMode>> // deduced execution mode
{};

TEST_P(ExecutionModeTest, call_update)
{
  auto exec_mode = tribol::getExecutionMode(std::get<0>(GetParam()), std::get<1>(GetParam()));
  EXPECT_EQ(exec_mode, std::get<2>(GetParam()));

  MPI_Barrier(MPI_COMM_WORLD);
}

INSTANTIATE_TEST_SUITE_P(tribol, ExecutionModeTest, testing::Values(
  std::make_tuple(tribol::MemorySpace::Host, tribol::ExecutionMode::Sequential, tribol::ExecutionMode::Sequential)
#ifdef TRIBOL_USE_OPENMP
  , std::make_tuple(tribol::MemorySpace::Host, tribol::ExecutionMode::Dynamic, tribol::ExecutionMode::OpenMP)
#else
  , std::make_tuple(tribol::MemorySpace::Host, tribol::ExecutionMode::Dynamic, tribol::ExecutionMode::Sequential)
#endif
#ifdef TRIBOL_USE_CUDA
  , std::make_tuple(tribol::MemorySpace::Host, tribol::ExecutionMode::Cuda, tribol::ExecutionMode::Sequential)
#endif
#ifdef TRIBOL_USE_HIP
  , std::make_tuple(tribol::MemorySpace::Host, tribol::ExecutionMode::Hip, tribol::ExecutionMode::Sequential)
#endif
#ifdef TRIBOL_USE_OPENMP
  , std::make_tuple(tribol::MemorySpace::Host, tribol::ExecutionMode::OpenMP, tribol::ExecutionMode::OpenMP)
#endif
  , std::make_tuple(tribol::MemorySpace::Dynamic, tribol::ExecutionMode::Sequential, tribol::ExecutionMode::Sequential)
  // should also throw warning
  , std::make_tuple(tribol::MemorySpace::Dynamic, tribol::ExecutionMode::Dynamic, tribol::ExecutionMode::Sequential)
#ifdef TRIBOL_USE_CUDA
  , std::make_tuple(tribol::MemorySpace::Dynamic, tribol::ExecutionMode::Cuda, tribol::ExecutionMode::Cuda)
#endif
#ifdef TRIBOL_USE_HIP
  , std::make_tuple(tribol::MemorySpace::Dynamic, tribol::ExecutionMode::Hip, tribol::ExecutionMode::Hip)
#endif
#ifdef TRIBOL_USE_OPENMP
  , std::make_tuple(tribol::MemorySpace::Dynamic, tribol::ExecutionMode::OpenMP, tribol::ExecutionMode::OpenMP)
#endif
#ifdef TRIBOL_USE_UMPIRE
  , std::make_tuple(tribol::MemorySpace::Device, tribol::ExecutionMode::Sequential, tribol::ExecutionMode::Sequential)
#if defined(TRIBOL_USE_CUDA)
  , std::make_tuple(tribol::MemorySpace::Device, tribol::ExecutionMode::Dynamic, tribol::ExecutionMode::Cuda)
#elif defined(TRIBOL_USE_HIP)
  , std::make_tuple(tribol::MemorySpace::Device, tribol::ExecutionMode::Dynamic, tribol::ExecutionMode::Hip)
#else
  // should also throw warning
  , std::make_tuple(tribol::MemorySpace::Device, tribol::ExecutionMode::Dynamic, tribol::ExecutionMode::Sequential)
#endif
#ifdef TRIBOL_USE_CUDA
  , std::make_tuple(tribol::MemorySpace::Device, tribol::ExecutionMode::Cuda, tribol::ExecutionMode::Cuda)
#endif
#ifdef TRIBOL_USE_HIP
  , std::make_tuple(tribol::MemorySpace::Device, tribol::ExecutionMode::Hip, tribol::ExecutionMode::Hip)
#endif
#ifdef TRIBOL_USE_OPENMP
  // should also throw warning
  , std::make_tuple(tribol::MemorySpace::Device, tribol::ExecutionMode::OpenMP, tribol::ExecutionMode::Sequential)
#endif
  , std::make_tuple(tribol::MemorySpace::Unified, tribol::ExecutionMode::Sequential, tribol::ExecutionMode::Sequential)
#if defined(TRIBOL_USE_CUDA)
  , std::make_tuple(tribol::MemorySpace::Unified, tribol::ExecutionMode::Dynamic, tribol::ExecutionMode::Cuda)
#elif defined(TRIBOL_USE_HIP)
  , std::make_tuple(tribol::MemorySpace::Unified, tribol::ExecutionMode::Dynamic, tribol::ExecutionMode::Hip)
#elif defined(TRIBOL_USE_OPENMP)
  , std::make_tuple(tribol::MemorySpace::Unified, tribol::ExecutionMode::Dynamic, tribol::ExecutionMode::OpenMP)
#else
  , std::make_tuple(tribol::MemorySpace::Unified, tribol::ExecutionMode::Dynamic, tribol::ExecutionMode::Sequential)
#endif
#ifdef TRIBOL_USE_CUDA
  , std::make_tuple(tribol::MemorySpace::Unified, tribol::ExecutionMode::Cuda, tribol::ExecutionMode::Cuda)
#endif
#ifdef TRIBOL_USE_HIP
  , std::make_tuple(tribol::MemorySpace::Unified, tribol::ExecutionMode::Hip, tribol::ExecutionMode::Hip)
#endif
#ifdef TRIBOL_USE_OPENMP
  , std::make_tuple(tribol::MemorySpace::Unified, tribol::ExecutionMode::OpenMP, tribol::ExecutionMode::OpenMP)
#endif
#endif
));

//------------------------------------------------------------------------------
#include "axom/slic/core/SimpleLogger.hpp"

int main(int argc, char* argv[])
{
  int result = 0;

  MPI_Init(&argc, &argv);

  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger;  // create & initialize test logger, finalized when
                                    // exiting main scope

  result = RUN_ALL_TESTS();

  MPI_Finalize();

  return result;
}
