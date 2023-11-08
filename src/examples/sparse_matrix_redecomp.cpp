// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Tribol Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (MIT)

/**
 * @file sparse_matrix_redecomp.cpp
 *
 * @brief Demonstrates moving a matrix using redecomp::SparseMatrixTransfer
 *
 * Demonstrates use of redecomp::SparseMatrixTransfer.
 *
 * For this example, the test space and trial space are the same, but
 * redecomp::SparseMatrixTransfer supports rectangular matrices, with different trial
 * and test spaces.
 *
 * Example runs (from repo root directory):
 *   - mpirun -np 4 {build_dir}/examples/sparse_matrix_redecomp_ex -r 1 -m
 *     data/two_hex.mesh
 *   - mpirun -np 4 {build_dir}/examples/sparse_matrix_redecomp_ex -r 1 -m
 *     data/star.mesh
 *   - mpirun -np 4 {build_dir}/examples/sparse_matrix_redecomp_ex -r 1 -o 2 -m
 *     data/two_hex.mesh
 *   - mpirun -np 4 {build_dir}/examples/sparse_matrix_redecomp_ex -r 1 -o 2 -m
 *     data/star.mesh
 */

#include <mfem/fem/fe_coll.hpp>
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

  double filter_radius = 1.0;

  axom::CLI::App app { "sparse_matrix_redecomp" };
  app.add_option("-m,--mesh", mesh_file, "Mesh file to use.")
    ->check(axom::CLI::ExistingFile)
    ->capture_default_str();
  app.add_option("-r,--refine", ref_levels,
    "Number of times to refine the mesh uniformly.")
    ->capture_default_str();
  app.add_option("-o,--order", order, 
    "Finite element order (polynomial degree).")
    ->capture_default_str();
  app.add_option("-R,--filter-radius", filter_radius,
    "Explicit filter radius.");
  CLI11_PARSE(app, argc, argv);

  SLIC_INFO_ROOT("Running sparse_matrix_redecomp with the following options:");
  SLIC_INFO_ROOT(axom::fmt::format("mesh:   {0}", mesh_file));
  SLIC_INFO_ROOT(axom::fmt::format("refine: {0}", ref_levels));
  SLIC_INFO_ROOT(axom::fmt::format("order:  {0}\n", order));

  SLIC_INFO_ROOT("Creating mfem::ParMesh...");
  // read serial mesh
  //auto mesh = std::make_unique<mfem::Mesh>(mesh_file.c_str(), 1, 1);
  auto mesh = mfem::Mesh::MakeCartesian2D(4, 4, mfem::Element::Type::QUADRILATERAL);

  // refine serial mesh
  for (int i{0}; i < ref_levels; ++i)
  {
    mesh.UniformRefinement();
  }
  
  // create parallel mesh from serial
  auto pmesh = std::make_unique<mfem::ParMesh>(MPI_COMM_WORLD, mesh);

  // further refinement of parallel mesh
  {
    int par_ref_levels = 0;
    for (int i{0}; i < par_ref_levels; ++i)
    {
      pmesh->UniformRefinement();
    }
  }

  mfem::VisItDataCollection pmesh_dc { "pmesh", pmesh.get() };

  // create finite element space on parmesh (1 point per element)
  mfem::L2_FECollection fe_coll { 0, pmesh->SpaceDimension() };
  mfem::ParFiniteElementSpace pmesh_space { pmesh.get(), &fe_coll };

  SLIC_INFO_ROOT("Creating redecomp::RedecompMesh...");
  // create redecomp mesh (ensure ghost region is >= filter radius)
  redecomp::RedecompMesh redecomp_mesh { *pmesh, redecomp::RedecompMesh::RCB, filter_radius };
  
  mfem::VisItDataCollection redecomp_dc { "redecomp_rank" + std::to_string(rank), &redecomp_mesh };
  
  // create finite element space on redecomp mesh (1 point per element)
  mfem::FiniteElementSpace redecomp_space { &redecomp_mesh, &fe_coll };

  // create transfer operator
  redecomp::SparseMatrixTransfer matrix_xfer {
    pmesh_space,
    pmesh_space,
    redecomp_space, 
    redecomp_space
  };

  SLIC_INFO_ROOT("Compute something on redecomp::RedecompMesh...");
  mfem::SparseMatrix smat { redecomp_space.GetVSize(), redecomp_space.GetVSize() };

  auto filter_kernel = [&filter_radius](const mfem::Vector &xi, const mfem::Vector &xj) {
    return std::max(0.0, filter_radius - xi.DistanceTo(xj));
  };

  int n_local_elem = redecomp_mesh.GetNE();
  int dim = redecomp_mesh.SpaceDimension();
  mfem::Vector xi(dim), xj(dim);
  std::vector<int> n_row_entries(n_local_elem);

  // loop over each element to form a row of the matrix
  for (int i=0; i<n_local_elem; i++) {
    redecomp_mesh.GetElementCenter(i, xi);

    // loop over each element to fill in each column (j) of row (i)
    n_row_entries[i] = 0;
    for (int j=0; j<n_local_elem; j++) {
      redecomp_mesh.GetElementCenter(j, xj);
      double Wij = filter_kernel(xi, xj);
      if (Wij > 0.0) {
        smat.Add(i, j, Wij);
        n_row_entries[i]++;
      }
    }
  }
  smat.Finalize();

  // normalize each row to conserve mass
  for (int i=0; i<n_local_elem; i++) {
    mfem::Vector row_data(smat.GetRowEntries(i), n_row_entries[i]);
    row_data /= row_data.Norml1();
  }

  SLIC_INFO_ROOT("Transferring something from RedecompMesh to ParMesh...");
  auto Mmat_xfer = matrix_xfer.TransferToParallel(smat);

  pmesh_dc.Save();
  redecomp_dc.Save();

  // test filter
  mfem::ParGridFunction x(&pmesh_space);
  mfem::ParGridFunction xf(&pmesh_space);

  mfem::FunctionCoefficient x_func([](const mfem::Vector&x) {
    return ((x[0] > 0.5) && (x[1] > 0.5))? 1.0 : 0.0;
  });
  
  x.ProjectCoefficient(x_func);
  
  Mmat_xfer->Mult(x, xf);

  mfem::VisItDataCollection visit_dc("test_filter", pmesh.get());
  visit_dc.RegisterField("x", &x);
  visit_dc.RegisterField("xf", &xf);
  visit_dc.Save();

  // cleanup
  MPI_Finalize();

  return 0;
}
