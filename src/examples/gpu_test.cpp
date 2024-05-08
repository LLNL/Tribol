// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Tribol Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (MIT)

#include <cstdio>

#include "examples_common.hpp" // for common functionality used in examples

#include "tribol/interface/tribol.hpp"

namespace tribol
{

struct ExampleMesh
{
  ArrayT<RealT, 2, MemorySpace::Device> coords;
  ArrayT<IndexT, 2, MemorySpace::Device> connectivity;
  ArrayT<RealT, 2, MemorySpace::Device> response;
};

TRIBOL_HOST_DEVICE void FillMesh(int i, RealT* coords, IndexT* conn, RealT* force)
{
  coords[i] = static_cast<RealT>(i);
  conn[i] = static_cast<IndexT>(i);
  force[i] = static_cast<RealT>(7 - i);
}

TRIBOL_HOST_DEVICE void OutputMesh(int i, RealT* coords)
{
  printf("Coordinate %d: %10.6f\n", i, coords[i]);
}

}

int main( int argc, char** argv )
{
  ////////////////////////////////
  //                            //
  // SETUP THE EXAMPLE AND MESH //
  //                            //
  ////////////////////////////////

  // initialize
#ifdef TRIBOL_USE_MPI
  MPI_Init( &argc, &argv );
#endif
  tribol::CommT problem_comm = TRIBOL_COMM_WORLD;
  initialize_logger( problem_comm );

  tribol::ExampleMesh mesh{{2, 4}, {2, 4}, {2, 4}};

  tribol::RealT* mesh_data = mesh.coords.data();
  tribol::IndexT* conn_data = mesh.connectivity.data();
  tribol::RealT* response_data = mesh.response.data();
  tribol::forAllExec<tribol::ExecutionMode::Cuda>(8, 
    [=] TRIBOL_HOST_DEVICE (tribol::IndexT i) {
      tribol::FillMesh(i, mesh_data, conn_data, response_data);
    }
  );
  tribol::forAllExec<tribol::ExecutionMode::Cuda>(8, 
    [=] TRIBOL_HOST_DEVICE (tribol::IndexT i) {
      tribol::OutputMesh(i, mesh_data);
    }
  );

  tribol::ArrayT<tribol::RealT, 2, tribol::MemorySpace::Host> host_coords(2, 4);
  host_coords(0, 0) = 2.0;
  host_coords(1, 0) = 3.0;
  host_coords(0, 1) = 4.0;
  host_coords(1, 1) = 5.0;
  host_coords(0, 2) = 6.0;
  host_coords(1, 2) = 7.0;
  host_coords(0, 3) = 8.0;
  host_coords(1, 3) = 9.0;
  mesh.coords = host_coords;
  mesh_data = mesh.coords.data();
  tribol::forAllExec<tribol::ExecutionMode::Cuda>(8, 
    [=] TRIBOL_HOST_DEVICE (int i) {
      tribol::OutputMesh(i, mesh_data);
    }
  );

  axom::slic::flushStreams();
  finalize_logger();

#ifdef TRIBOL_USE_MPI
  MPI_Finalize();
#endif

  return 0;
}

