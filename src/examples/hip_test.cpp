// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Tribol Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (MIT)

#include <iostream>

#include "tribol/common/ArrayTypes.hpp"
#include "tribol/common/ExecModel.hpp"


int main(int argc, char** argv)
{
  auto allocator_id = axom::getUmpireResourceAllocatorID(umpire::resource::MemoryResourceType::Device);
  axom::Array<int, 1, axom::MemorySpace::Dynamic> test({71}, allocator_id);

  allocator_id = axom::getUmpireResourceAllocatorID(umpire::resource::MemoryResourceType::Host);
  axom::Array<int, 1, axom::MemorySpace::Dynamic> test_host(test, allocator_id);

  std::cout << test_host[0] << std::endl;

  axom::Array<int, 1, axom::MemorySpace::Host> test_host2(test);

  std::cout << test_host2[0] << std::endl;

  return 0;
}
