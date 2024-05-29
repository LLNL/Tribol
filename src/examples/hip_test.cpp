// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Tribol Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (MIT)

#include <iostream>
#include <cstdio>

#include "tribol/common/ArrayTypes.hpp"
#include "tribol/common/ExecModel.hpp"
#include "tribol/common/LoopExec.hpp"

int main(int argc, char** argv)
{
  auto allocator_id = axom::getUmpireResourceAllocatorID(umpire::resource::MemoryResourceType::Device);
  axom::Array<int> test({71}, allocator_id);

  //hipDeviceSynchronize();

  std::cout << "On host, before copy:" << test[0] << std::endl;

  auto test_view = test.view();
  std::cout << "On host, after view:" << test[0] << std::endl;
  forAllExec(tribol::ExecutionMode::Hip, 1, [test_view] TRIBOL_HOST_DEVICE (int i) mutable {
    printf("On device: %d\n", test_view[i]);
  });
  std::cout << "On host, after print device:" << test[0] << std::endl;


  allocator_id = axom::getUmpireResourceAllocatorID(umpire::resource::MemoryResourceType::Host);
  axom::Array<int, 1, axom::MemorySpace::Dynamic> test_host(test, allocator_id);

  std::cout << "On host, after copy:" << test[0] << std::endl;

  //hipDeviceSynchronize();
  std::cout << test_host[0] << std::endl;

  axom::Array<int, 1, axom::MemorySpace::Host> test_host2(test);

  //hipDeviceSynchronize();

  std::cout << test_host2[0] << std::endl;

  // create a device array with some random data
  allocator_id = axom::getUmpireResourceAllocatorID(umpire::resource::MemoryResourceType::Device);
  axom::Array<double> test2({12.0, 21.0, 101.0, 110.0}, allocator_id);

  // hipDeviceSynchronize();

  // create a view of that data held in an array on host
  axom::Array<axom::ArrayView<double>> tests_host(1, 1);
  tests_host[0] = axom::ArrayView<double>(test2.data(), 4);

  // hipDeviceSynchronize();

  // copy that array of the view to the device
  axom::Array<axom::ArrayView<double>> tests(tests_host, allocator_id);

  // hipDeviceSynchronize();

  // create a view of that array of the view
  axom::ArrayView<axom::ArrayView<double>> tests_view = tests.view();

  // hipDeviceSynchronize();

  // print them on device
  forAllExec(tribol::ExecutionMode::Hip, 4, [tests_view] TRIBOL_HOST_DEVICE (int i) {
    printf("tests[0][%d]: %10.4f\n", i, tests_view[0][i]);
  });

  axom::Array<double, 2> test3({3, 1}, allocator_id);
  test3.fill(0.0);

  // print them on device
  auto test3_view = test3.view();
  forAllExec(tribol::ExecutionMode::Hip, 3, [test3_view] TRIBOL_HOST_DEVICE (int i) {
    printf("test3(%d,0): %10.4f\n", i, test3_view(i,0));
  });

  return 0;
}
