// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Tribol Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (MIT)

/**
 * @file matrix_redecomp.cpp
 *
 * @brief Demonstrates moving a matrix using redecomp::MatrixTransfer
 *
 * Demonstrates use of redecomp::MatrixTransfer to move element-level matrices
 * from finite element spaces on a redecomp::RedecompMesh to an assembled,
 * global matrix on finite element spaces on a mfem::ParMesh. Specifically,
 * element contributions of the mass matrix are computed on the
 * redecomp::RedecompMesh, then transferred to the parent mfem::ParMesh where
 * they are assembled.  These mass matrix values are compared to mass matrix
 * values computed directly on the mfem::ParMesh.
 *
 * For this example, the test space and trial space are the same, but
 * redecomp::MatrixTransfer supports rectangular matrices, with different trial
 * and test spaces.
 *
 * Example runs (from repo root directory):
 *   - mpirun -np 4 {build_dir}/examples/matrix_redecomp_ex -r 1 -m
 *     data/two_hex.mesh
 *   - mpirun -np 4 {build_dir}/examples/matrix_redecomp_ex -r 1 -m
 *     data/star.mesh
 *   - mpirun -np 4 {build_dir}/examples/matrix_redecomp_ex -r 1 -o 2 -m
 *     data/two_hex.mesh
 *   - mpirun -np 4 {build_dir}/examples/matrix_redecomp_ex -r 1 -o 2 -m
 *     data/star.mesh
 */

#include <string>

#include <mpi.h>

#include "axom/CLI11.hpp"
#include "axom/slic.hpp"
#include "mfem.hpp"

#include "redecomp/redecomp.hpp"
#include "redecomp/utils/ArrayUtility.hpp"
#include "tribol/config.hpp"

int main( int argc, char** argv )
{
  // initialize MPI
  MPI_Init( &argc, &argv );
  int np, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // initialize logger
  axom::slic::SimpleLogger logger;
  axom::slic::setIsRoot(rank == 0);

  // command line options
  // location of mesh file. TRIBOL_REPO_DIR is defined in tribol/config.hpp
  std::string mesh_file = TRIBOL_REPO_DIR "/data/star.mesh";
  // number of times to uniformly refine the serial mesh before constructing the
  // parallel mesh
  int ref_levels = 0;
  // polynomial order of the finite element discretization
  int order = 1;

  axom::CLI::App app { "matrix_redecomp" };
  app.add_option("-m,--mesh", mesh_file, "Mesh file to use.")
    ->check(axom::CLI::ExistingFile)
    ->capture_default_str();
  app.add_option("-r,--refine", ref_levels,
    "Number of times to refine the mesh uniformly.")
    ->capture_default_str();
  app.add_option("-o,--order", order, 
    "Finite element order (polynomial degree).")
    ->capture_default_str();
  CLI11_PARSE(app, argc, argv);

  SLIC_INFO_ROOT("Running quadrature_matrix_redecomp with the following options:");
  SLIC_INFO_ROOT(axom::fmt::format("mesh:   {0}", mesh_file));
  SLIC_INFO_ROOT(axom::fmt::format("refine: {0}", ref_levels));
  SLIC_INFO_ROOT(axom::fmt::format("order:  {0}\n", order));

  SLIC_INFO_ROOT("Creating mfem::ParMesh...");
  // read serial mesh
  auto mesh = std::make_unique<mfem::Mesh>(mesh_file.c_str(), 1, 1);

  // refine serial mesh
  for (int i{0}; i < ref_levels; ++i)
  {
    mesh->UniformRefinement();
  }
  
  // create parallel mesh from serial
  auto pmesh = std::make_unique<mfem::ParMesh>(MPI_COMM_WORLD, *mesh);
  mesh.reset(nullptr);

  // further refinement of parallel mesh
  {
    int par_ref_levels = 0;
    for (int i{0}; i < par_ref_levels; ++i)
    {
      pmesh->UniformRefinement();
    }
  }

  mfem::VisItDataCollection pmesh_dc { "pmesh", pmesh.get() };

  // create quadrature space on parmesh (1 point per element)
  mfem::QuadratureSpace pmesh_quad_space { pmesh.get(), 1 };

  SLIC_INFO_ROOT("Creating redecomp::RedecompMesh...");
  // create redecomp mesh
  redecomp::RedecompMesh redecomp_mesh { *pmesh, redecomp::RedecompMesh::RCB, 1.0 };
  
  mfem::VisItDataCollection redecomp_dc { "redecomp", &redecomp_mesh };
  
  // create quadrature space on redecomp mesh (1 point per element)
  mfem::QuadratureSpace redecomp_quad_space { &redecomp_mesh, 1 };

  // create transfer operator
  redecomp::QuadratureMatrixTransfer matrix_xfer {
    pmesh_quad_space,
    pmesh_quad_space,
    redecomp_quad_space, 
    redecomp_quad_space
  };

  SLIC_INFO_ROOT("Compute something on redecomp::RedecompMesh...");
  mfem::SparseMatrix smat { redecomp_quad_space.GetSize(), redecomp_quad_space.GetSize() };
  // Compute something here...
  smat.Finalize();

  SLIC_INFO_ROOT("Transferring something from RedecompMesh to ParMesh...");
  auto Mmat_xfer = matrix_xfer.TransferToParallel(smat);

  pmesh_dc.Save();
  redecomp_dc.Save();

  // cleanup
  MPI_Finalize();

  return 0;
}
