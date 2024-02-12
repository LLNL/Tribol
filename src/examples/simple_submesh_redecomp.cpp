// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Tribol Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (MIT)

/**
 * @file simple_submesh_redecomp.cpp
 *
 * @brief Demonstrates redecomp::RedecompMesh and redecomp transfer capability
 *
 * Demonstrates use of redecomp::RedecompMesh to build a RCB-based
 * redecomposition of an mfem::ParMesh, then transfers an mfem::ParGridFunction
 * and an mfem::QuadratureFunction to an mfem::GridFunction and an
 * mfem::QuadratureFunction on the redecomp::RedecompMesh and vice versa.
 *
 * Though redecomp::RedecompMesh contains pieces of the mesh spread across
 * ranks, it is a serial mesh which derives from mfem::Mesh. All coordination
 * with other ranks must be done by transferring data to the parent
 * mfem::ParMesh level. Beyond the decomposition method (RCB vs. typically metis
 * k-way), a redecomp::RedecompMesh differs from an mfem::ParMesh because it
 * includes a layer of ghost elements at the on-rank domain boundaries.
 * Quantities on ghost elements are typically transferred to fields on the
 * redecomp::RedecompMesh, but are not transferred back to the mfem::ParMesh.
 *
 * This example illustrates the envisioned typical workflow for using
 * redecomp::RedecompMesh: 
 *   1. Create a RedecompMesh by providing a parent ParMesh (and optionally a
 *      method of repartitioning the domain)
 *   2. Create fields on the RedecompMesh (mfem::GridFunction and/or
 *      mfem::QuadratureFunction)
 *   3. Creating a RedecompTransfer object which provides methods to transfer an
 *      mfem::ParGridFunction and an mfem::QuadratureFunction from the parent
 *      ParMesh to the RedecompMesh and vice versa
 *   4. Transferring the mfem::ParGridfunction and/or mfem::QuadratureFunction
 *      on the parent ParMesh to the mfem::GridFunction and/or
 *      mfem::QuadratureFunction on the RedecompMesh
 *   5. Performing an update on the mfem::GridFunction and/or
 *      mfem::QuadratureFunction at the RedecompMesh level
 *   6. Transferring the mfem::GridFunction and/or mfem::QuadratureFunction on
 *      the RedecompMesh back to the corresponding functions on the ParMesh
 *
 * Example runs (from repo root directory):
 *   - mpirun -np 4 {build_dir}/examples/domain_redecomp_ex -r 1 -m data/two_hex.mesh
 *   - mpirun -np 4 {build_dir}/examples/domain_redecomp_ex -r 1 -m data/star.mesh
 *   - mpirun -np 4 {build_dir}/examples/domain_redecomp_ex -r 1 -o 2 -m data/two_hex.mesh
 *   - mpirun -np 4 {build_dir}/examples/domain_redecomp_ex -r 1 -o 2 -m data/star.mesh
 *
 * Example output on both the RedecompMesh and the ParMesh can be viewed in
 * VisIt.
 */

#include <mpi.h>

// MFEM includes
#include "mfem.hpp"

// Axom includes
#include "axom/CLI11.hpp"
#include "axom/slic.hpp"

// Redecomp includes
#include "redecomp/redecomp.hpp"

// Tribol includes
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
  // radius of the ghost layer of elements
  double ghost_radius = 0.05;

  axom::CLI::App app { "simple_submesh_redecomp" };
  app.add_option("-m,--mesh", mesh_file, "Mesh file to use.")
    ->check(axom::CLI::ExistingFile)
    ->capture_default_str();
  app.add_option("-r,--refine", ref_levels,
    "Number of times to refine the mesh uniformly.")
    ->capture_default_str();
  app.add_option("-o,--order", order, 
    "Finite element order (polynomial degree).")
    ->capture_default_str();
  app.add_option("-g,--ghost-radius", ghost_radius,
    "Radius of the ghost layer of elements.")
    ->capture_default_str();
  CLI11_PARSE(app, argc, argv);

  SLIC_INFO_ROOT("Running simple_submesh_redecomp with the following options:");
  SLIC_INFO_ROOT(axom::fmt::format("mesh:         {0}", mesh_file));
  SLIC_INFO_ROOT(axom::fmt::format("refine:       {0}", ref_levels));
  SLIC_INFO_ROOT(axom::fmt::format("order:        {0}", order));
  SLIC_INFO_ROOT(axom::fmt::format("ghost radius: {0}\n", ghost_radius));

  // Read the mesh, refine, and create a mfem::ParMesh
  mfem::Mesh serial_mesh(mesh_file);
  for (int i = 0; i < ref_levels; ++i)
  {
    serial_mesh.UniformRefinement();
  }
  mfem::ParMesh mesh(MPI_COMM_WORLD, serial_mesh);
  serial_mesh.Clear();
  // Further refinement of parallel mesh
  int par_ref_levels = 2;
  for (int i = 0; i < par_ref_levels; ++i)
  {
    mesh.UniformRefinement();
  }

  // Create an H1 finite element space on the mesh for coordinates
  mfem::H1_FECollection fec(order, mesh.Dimension());
  mfem::ParFiniteElementSpace fespace(&mesh, &fec, mesh.SpaceDimension());
  mfem::ParGridFunction coords(&fespace);
  mesh.SetNodalGridFunction(&coords);

  // Output mesh and grid function
  mfem::VisItDataCollection mesh_dc("mesh", &mesh);
  mesh_dc.RegisterField("coordinates", &coords);
  mesh_dc.Save();

  // Create a surface submesh
  mfem::ParSubMesh submesh = mfem::ParSubMesh::CreateFromBoundary(mesh, mesh.bdr_attributes);

  // Create a redecomp mesh
  redecomp::RedecompMesh redecomp_mesh(submesh, ghost_radius);

  // Move the coordinates grid function to submesh manually
  mfem::H1_FECollection submesh_fec(order, submesh.Dimension());
  mfem::ParFiniteElementSpace submesh_fespace(&submesh, &submesh_fec, submesh.SpaceDimension());
  mfem::ParGridFunction submesh_coords(&submesh_fespace);
  submesh.Transfer(coords, submesh_coords);

  // Output submesh and grid function
  mfem::VisItDataCollection submesh_dc("submesh", &submesh);
  submesh_dc.RegisterField("coordinates", &submesh_coords);
  submesh_dc.Save();

  // Move submesh coordinates grid function to redecomp mesh manually
  mfem::H1_FECollection redecomp_fec(order, redecomp_mesh.Dimension());
  mfem::FiniteElementSpace redecomp_fespace(&redecomp_mesh, &redecomp_fec, redecomp_mesh.SpaceDimension());
  mfem::GridFunction redecomp_coords(&redecomp_fespace);
  redecomp::RedecompTransfer redecomp_xfer;
  redecomp_xfer.TransferToSerial(submesh_coords, redecomp_coords);

  // Output redecomp mesh and grid function
  mfem::VisItDataCollection redecomp_dc("redecomp_" + std::to_string(rank), &redecomp_mesh);
  redecomp_dc.RegisterField("coordinates", &redecomp_coords);
  redecomp_dc.Save();

  // cleanup
  MPI_Finalize();

  return 0;
}
