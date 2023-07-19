// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Tribol Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (MIT)

#include <set>

#ifdef TRIBOL_USE_UMPIRE
// Umpire includes
#include "umpire/ResourceManager.hpp"
#endif

// MFEM includes
#include "mfem.hpp"

// Axom includes
#include "axom/CLI11.hpp"
#include "axom/slic.hpp"

// Redecomp includes
#include "redecomp/redecomp.hpp"

// Tribol includes
#include "tribol/common/Parameters.hpp"
#include "tribol/config.hpp"
#include "tribol/interface/tribol.hpp"
#include "tribol/interface/mfem_tribol.hpp"


/**
 * @file mfem_mortar_lm_patch.cpp
 *
 * @brief Demonstrates contact patch test using the mortar method
 *
 * Demonstrates a contact patch test using the mortar method in Tribol.
 * Enforcement is through Lagrange multipliers and no active set (i.e. tied +
 * sliding contact).  Small deformation contact is assumed and, consequently,
 * the solution is determined through a single linear solve (no timestepping).
 *
 * The example uses the Tribol MFEM interface, which supports decomposed (MPI)
 * meshes and will support higher order meshes (through LOR) in a future update.
 *
 * Example runs (from repo root directory):
 *   - mpirun -np 4 {build_dir}/examples/mfem_mortar_lm_patch_ex
 *
 * Example output can be viewed in VisIt or ParaView.
 */
int main( int argc, char** argv )
{
  // initialize MPI
  MPI_Init( &argc, &argv );
  int np, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

#ifdef TRIBOL_USE_UMPIRE
  umpire::ResourceManager::getInstance();  // initialize umpire's ResouceManager
#endif

  // initialize logger
  axom::slic::SimpleLogger logger;
  axom::slic::setIsRoot(rank == 0);

  // command line options
  // number of times to uniformly refine the serial mesh before constructing the
  // parallel mesh
  int ref_levels = 2;
  // polynomial order of the finite element discretization
  int order = 1;
  // Lame parameter
  double lambda = 50.0;
  // Lame parameter (shear modulus)
  double mu = 50.0;

  axom::CLI::App app { "mfem_mortar_lm_patch" };
  app.add_option("-r,--refine", ref_levels,
    "Number of times to refine the mesh uniformly.")
    ->capture_default_str();
  // TODO: LOR support for implicit contact
  // app.add_option("-o,--order", order, 
  //   "Finite element order (polynomial degree).")
  //   ->capture_default_str();
  app.add_option("-l,--lambda", lambda, 
    "Lame parameter.")
    ->capture_default_str();
  app.add_option("-m,--mu", mu, 
    "Lame parameter (shear modulus).")
    ->capture_default_str();
  CLI11_PARSE(app, argc, argv);

  SLIC_INFO_ROOT("Running mfem_mortar_lm_patch with the following options:");
  SLIC_INFO_ROOT(axom::fmt::format("refine: {0}", ref_levels));
  SLIC_INFO_ROOT(axom::fmt::format("order:  {0}\n", order));

  // fixed options
  // location of mesh file. TRIBOL_REPO_DIR is defined in tribol/config.hpp
  std::string mesh_file = TRIBOL_REPO_DIR "/data/two_hex_overlap.mesh";
  // boundary element attributes of mortar surface
  auto mortar_attrs = std::set<int>({4});
  // boundary element attributes of nonmortar surface 
  auto nonmortar_attrs = std::set<int>({5});
  // boundary element attributes of x-fixed surfaces (at x = 0)
  auto xfixed_attrs = std::set<int>({1});
  // boundary element attributes of y-fixed surfaces (at y = 0)
  auto yfixed_attrs = std::set<int>({2});
  // boundary element attributes of z-fixed surfaces (3: surface at z = 0, 6: surface at z = 1.95)
  auto zfix_attribs = std::set<int>({3, 6});
  
  axom::utilities::Timer timer { false };

  // read mesh
  timer.start();
  std::unique_ptr<mfem::ParMesh> pmesh { nullptr };
  {
    // read serial mesh
    auto mesh = std::make_unique<mfem::Mesh>(mesh_file.c_str(), 1, 1);

    // refine serial mesh
    if (ref_levels > 0)
    {
      for (int i{0}; i < ref_levels; ++i)
      {
        mesh->UniformRefinement();
      }
    }
    
    // create parallel mesh from serial
    pmesh = std::make_unique<mfem::ParMesh>(MPI_COMM_WORLD, *mesh);
    mesh.reset(nullptr);

    // further refinement of parallel mesh
    {
      int par_ref_levels = 0;
      for (int i{0}; i < par_ref_levels; ++i)
      {
        pmesh->UniformRefinement();
      }
    }
  }
  timer.stop();
  SLIC_INFO_ROOT(axom::fmt::format(
    "Time to create parallel mesh: {0:f}ms", timer.elapsedTimeInMilliSec()
  ));
  
  // set up data collection for output
  auto paraview_datacoll = mfem::ParaViewDataCollection("mortar_patch_pv", pmesh.get());
  auto visit_datacoll = mfem::VisItDataCollection("mortar_patch_vi", pmesh.get());

  // grid function for higher-order nodes
  timer.start();
  auto fe_coll = mfem::H1_FECollection(order, pmesh->SpaceDimension());
  auto par_fe_space = mfem::ParFiniteElementSpace(
    pmesh.get(), &fe_coll, pmesh->SpaceDimension());
  auto coords = mfem::ParGridFunction(&par_fe_space);
  if (order > 1)
  {
    pmesh->SetNodalGridFunction(&coords, false);
  }
  else
  {
    pmesh->GetNodes(coords);
  }
  paraview_datacoll.RegisterField("pos", &coords);
  visit_datacoll.RegisterField("pos", &coords);

  // grid function for displacement
  mfem::ParGridFunction displacement { &par_fe_space };
  paraview_datacoll.RegisterField("disp", &displacement);
  visit_datacoll.RegisterField("disp", &displacement);
  displacement = 0.0;
  timer.stop();
  SLIC_INFO_ROOT(axom::fmt::format(
    "Time to create grid functions: {0:f}ms", timer.elapsedTimeInMilliSec()
  ));

  // save initial configuration
  paraview_datacoll.Save();
  visit_datacoll.Save();

  // recover dirichlet bc tdof list
  timer.start();
  mfem::Array<int> ess_tdof_list;
  {
    mfem::Array<int> ess_vdof_marker;
    mfem::Array<int> ess_bdr(pmesh->bdr_attributes.Max());
    ess_bdr = 0;
    for (auto xfixed_attr : xfixed_attrs)
    {
      ess_bdr[xfixed_attr-1] = 1;
    }
    par_fe_space.GetEssentialVDofs(ess_bdr, ess_vdof_marker, 0);
    mfem::Array<int> new_ess_vdof_marker;
    ess_bdr = 0;
    for (auto yfixed_attr : yfixed_attrs)
    {
      ess_bdr[yfixed_attr-1] = 1;
    }
    par_fe_space.GetEssentialVDofs(ess_bdr, new_ess_vdof_marker, 1);
    for (int i{0}; i < ess_vdof_marker.Size(); ++i)
    {
      ess_vdof_marker[i] = ess_vdof_marker[i] || new_ess_vdof_marker[i];
    }
    ess_bdr = 0;
    for (auto zfix_attrib : zfix_attribs)
    {
      ess_bdr[zfix_attrib-1] = 1;
    }
    par_fe_space.GetEssentialVDofs(ess_bdr, new_ess_vdof_marker, 2);
    for (int i{0}; i < ess_vdof_marker.Size(); ++i)
    {
      ess_vdof_marker[i] = ess_vdof_marker[i] || new_ess_vdof_marker[i];
    }
    mfem::Array<int> ess_tdof_marker;
    par_fe_space.GetRestrictionMatrix()->BooleanMult(ess_vdof_marker, ess_tdof_marker);
    mfem::FiniteElementSpace::MarkerToList(ess_tdof_marker, ess_tdof_list);
  }
  timer.stop();
  SLIC_INFO_ROOT(axom::fmt::format(
    "Time to set up boundary conditions: {0:f}ms", timer.elapsedTimeInMilliSec()
  ));

  // set up mfem elasticity bilinear form
  timer.start();
  mfem::ParBilinearForm a(&par_fe_space);
  mfem::ConstantCoefficient lambda_coeff(lambda);
  mfem::ConstantCoefficient mu_coeff(mu);
  a.AddDomainIntegrator(new mfem::ElasticityIntegrator(lambda_coeff, mu_coeff));

  // compute elasticity contribution to stiffness
  a.Assemble();
  auto A = std::make_unique<mfem::HypreParMatrix>();
  a.FormSystemMatrix(ess_tdof_list, *A);
  timer.stop();
  SLIC_INFO_ROOT(axom::fmt::format(
    "Time to create and assemble internal stiffness: {0:f}ms", timer.elapsedTimeInMilliSec()
  ));

  // set up tribol
  timer.start();
  tribol::initialize(pmesh->SpaceDimension(), MPI_COMM_WORLD);
  int coupling_scheme_id = 0;
  int mesh1_id = 0;
  int mesh2_id = 1;
  tribol::registerMfemCouplingScheme(
    coupling_scheme_id, mesh1_id, mesh2_id,
    *pmesh, coords, mortar_attrs, nonmortar_attrs,
    tribol::SURFACE_TO_SURFACE, 
    tribol::NO_SLIDING, 
    tribol::SINGLE_MORTAR, 
    tribol::FRICTIONLESS,
    tribol::LAGRANGE_MULTIPLIER,
    tribol::BINNING_GRID
  );
  tribol::setLagrangeMultiplierOptions(
    0,
    tribol::ImplicitEvalMode::MORTAR_RESIDUAL_JACOBIAN
  );

  // update tribol (compute contact contribution to force and stiffness)
  tribol::updateMfemParallelDecomposition();
  double dt {1.0};  // time is arbitrary here (no timesteps)
  tribol::update(1, 1.0, dt);
  paraview_datacoll.SetCycle(1);
  paraview_datacoll.SetTime(1.0);
  paraview_datacoll.SetTimeStep(1.0);
  visit_datacoll.SetCycle(1);
  visit_datacoll.SetTime(1.0);
  visit_datacoll.SetTimeStep(1.0);

  // retrieve block stiffness matrix
  auto A_blk = tribol::getMfemBlockJacobian(0);
  A_blk->SetBlock(0, 0, A.release());
  timer.stop();
  SLIC_INFO_ROOT(axom::fmt::format(
    "Time to setup Tribol and compute Jacobian: {0:f}ms", timer.elapsedTimeInMilliSec()
  ));

  // create block solution and RHS vectors
  timer.start();
  mfem::BlockVector B_blk { A_blk->ColOffsets() };
  B_blk = 0.0;
  mfem::BlockVector X_blk { A_blk->RowOffsets() };
  X_blk = 0.0;

  // retrieve gap vector (RHS) from contact
  mfem::ParGridFunction g;
  tribol::getMfemGap(0, g);

  // restriction operator on submesh: maps dofs stored in g to tdofs stored in G
  {
    auto& G = B_blk.GetBlock(1);
    auto& R_submesh = *g.ParFESpace()->GetRestrictionOperator();
    R_submesh.Mult(g, G);
  }

  // solve for X_blk
  mfem::MINRESSolver solver(MPI_COMM_WORLD);
  solver.SetRelTol(1.0e-8);
  solver.SetAbsTol(1.0e-12);
  solver.SetMaxIter(5000);
  solver.SetPrintLevel(3);
  solver.SetOperator(*A_blk);
  solver.Mult(B_blk, X_blk);

  // move block displacements to grid function
  {
    auto& U = X_blk.GetBlock(0);
    auto& P = *par_fe_space.GetProlongationMatrix();
    P.Mult(U, displacement);
  }
  displacement.Neg();

  // update mesh coordinates
  coords += displacement;
  pmesh->SetVertices(coords);
  timer.stop();
  SLIC_INFO_ROOT(axom::fmt::format(
    "Time to solve for updated displacements: {0:f}ms", timer.elapsedTimeInMilliSec()
  ));

  // save deformed configuration
  paraview_datacoll.Save();
  visit_datacoll.Save();

  // cleanup
  tribol::finalize();
  MPI_Finalize();

  return 0;
}
