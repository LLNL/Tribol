// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Tribol Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (MIT)

/**
 * @file mfem_mortar_lm_patch.cpp
 *
 * @brief Demonstrates contact patch test using the mortar method
 *
 * Demonstrates a three dimensional contact patch test using the mortar method
 * in Tribol. Contact is enforced between two blocks which are initially in
 * contact. The blocks occupy [0, 1]^3 and [0, 1]x[0, 1]x[0.99, 1.99]. To
 * enforce symmetry and prevent rigid body modes, dirichlet boundary conditions
 * are applied in the x-direction along the x = 0 plane, in the y-direction
 * along y = 0 plane, and in the z-direction along the z = 0 and z = 1.99
 * planes. Enforcement is through Lagrange multipliers and no active set (i.e.
 * tied + sliding contact). Small deformation contact is assumed and,
 * consequently, the system is linear and the solution is determined through a
 * single linear solve (no timestepping).
 *
 * The linear system solved is
 *  | A B^T | | d | = | f |
 *  | B 0   | | p |   | g | ,
 * where A is the system matrix for elasticity, B is the constraint matrix for
 * mortar contact, d is the vector of nodal displacements, p is the vector of
 * contact pressures, f is the vector of nodal forces, and g is the vector of
 * gaps on the contact surface.  MFEM block operators and vectors are used to 
 * store the linear system.
 *
 * The example uses the Tribol MFEM interface, which supports decomposed (MPI)
 * meshes and will support higher order meshes (through LOR) in a future update.
 * Comments in the main function below give details on each step of the example
 * code.
 *
 * Example runs (from repo root directory):
 *   - mpirun -np 4 {build_dir}/examples/mfem_mortar_lm_patch_ex
 *
 * Example output can be viewed in VisIt or ParaView.
 */

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

  // define command line options
  // number of times to uniformly refine the serial mesh before constructing the
  // parallel mesh
  int ref_levels = 2;
  // polynomial order of the finite element discretization
  int order = 1;
  // Lame parameter lambda
  double lambda = 50.0;
  // Lame parameter mu (shear modulus)
  double mu = 50.0;

  // parse command line options
  axom::CLI::App app { "mfem_mortar_lm_patch" };
  app.add_option("-r,--refine", ref_levels,
    "Number of times to refine the mesh uniformly.")
    ->capture_default_str();
  // TODO: LOR support for implicit contact
  // app.add_option("-o,--order", order, 
  //   "Finite element order (polynomial degree).")
  //   ->capture_default_str();
  app.add_option("-l,--lambda", lambda, 
    "Lame parameter lambda.")
    ->capture_default_str();
  app.add_option("-m,--mu", mu, 
    "Lame parameter mu (shear modulus).")
    ->capture_default_str();
  CLI11_PARSE(app, argc, argv);

  SLIC_INFO_ROOT("Running mfem_mortar_lm_patch with the following options:");
  SLIC_INFO_ROOT(axom::fmt::format("refine: {0}", ref_levels));
  SLIC_INFO_ROOT(axom::fmt::format("lambda: {0}", lambda));
  SLIC_INFO_ROOT(axom::fmt::format("mu:     {0}\n", mu));

  // fixed options
  // location of mesh file. TRIBOL_REPO_DIR is defined in tribol/config.hpp
  std::string mesh_file = TRIBOL_REPO_DIR "/data/two_hex_overlap.mesh";
  // boundary element attributes of mortar surface, the z = 1 plane of the first
  // block
  auto mortar_attrs = std::set<int>({4});
  // boundary element attributes of nonmortar surface, the z = 0.99 plane of the
  // second block
  auto nonmortar_attrs = std::set<int>({5});
  // boundary element attributes of x-fixed surfaces (at x = 0)
  auto xfixed_attrs = std::set<int>({1});
  // boundary element attributes of y-fixed surfaces (at y = 0)
  auto yfixed_attrs = std::set<int>({2});
  // boundary element attributes of z-fixed surfaces (3: surface at z = 0, 6:
  // surface at z = 1.99)
  auto zfix_attribs = std::set<int>({3, 6});
  
  // create an axom timer to give wall times for each step
  axom::utilities::Timer timer { false };

  // This block of code will read the mesh data given in two_hex_overlap.mesh,
  // create an mfem::Mesh, refine the mesh, then create an mfem::ParMesh.
  // Optionally, the mfem::ParMesh can be refined further on each rank by
  // setting par_ref_levels >= 1, though this is disabled below.
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
      // set this to >= 1 to refine the mesh on each rank further
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
  
  // Set up an MFEM data collection for output. We output data in Paraview and
  // VisIt formats.
  auto paraview_datacoll = mfem::ParaViewDataCollection("mortar_patch_pv", pmesh.get());
  auto visit_datacoll = mfem::VisItDataCollection("mortar_patch_vi", pmesh.get());

  // This block of code creates position and displacement grid functions (and
  // associated finite element collections and finite element spaces) on the
  // mesh. The displacement grid function is initialized to zero. The position
  // grid function is initialized with the nodal coordinates. These grid
  // functions are registered with the data collections for output.
  timer.start();
  // Finite element collection (shared between all grid functions).
  auto fe_coll = mfem::H1_FECollection(order, pmesh->SpaceDimension());
  // Finite element space (shared between all grid functions).
  auto par_fe_space = mfem::ParFiniteElementSpace(
    pmesh.get(), &fe_coll, pmesh->SpaceDimension());
  // Create coordinate grid function
  auto coords = mfem::ParGridFunction(&par_fe_space);
  // Set coordinate grid function based on nodal locations. In MFEM, nodal
  // locations of higher order meshes are stored in a grid function. For linear
  // MFEM meshes, nodal locations can be stored in a grid function or through
  // the vertex coordinates. For consistency, we will create a nodal grid
  // function even for linear meshes.
  pmesh->SetNodalGridFunction(&coords, false);
  paraview_datacoll.RegisterField("pos", &coords);
  visit_datacoll.RegisterField("pos", &coords);

  // Create a grid function for displacement
  mfem::ParGridFunction displacement { &par_fe_space };
  paraview_datacoll.RegisterField("disp", &displacement);
  visit_datacoll.RegisterField("disp", &displacement);
  displacement = 0.0;
  timer.stop();
  SLIC_INFO_ROOT(axom::fmt::format(
    "Time to create grid functions: {0:f}ms", timer.elapsedTimeInMilliSec()
  ));

  // This block of code builds a list of degrees of freedom to which homogeneous
  // displacement boundary conditions will be applied. These boundary conditions
  // enforce problem symmetry and prevent rigid body deformation modes. The
  // boundary attribute sets for each direction are identified in the fixed
  // options above.
  timer.start();
  mfem::Array<int> ess_tdof_list;
  {
    // First, build an array of "markers" (i.e. booleans) to denote which vdofs
    // are in the list.
    mfem::Array<int> ess_vdof_marker;
    // Convert x-fixed boundary attributes into markers
    mfem::Array<int> ess_bdr(pmesh->bdr_attributes.Max());
    ess_bdr = 0;
    for (auto xfixed_attr : xfixed_attrs)
    {
      ess_bdr[xfixed_attr-1] = 1;
    }
    // Find all x-direction vdofs with the given boundary marker
    par_fe_space.GetEssentialVDofs(ess_bdr, ess_vdof_marker, 0);
    // Convert y-fixed boundary attributes into markers
    mfem::Array<int> new_ess_vdof_marker;
    ess_bdr = 0;
    for (auto yfixed_attr : yfixed_attrs)
    {
      ess_bdr[yfixed_attr-1] = 1;
    }
    // Find all y-direction vdofs with the given boundary marker
    par_fe_space.GetEssentialVDofs(ess_bdr, new_ess_vdof_marker, 1);
    // Compute union of x-direction and y-direction vdofs
    for (int i{0}; i < ess_vdof_marker.Size(); ++i)
    {
      ess_vdof_marker[i] = ess_vdof_marker[i] || new_ess_vdof_marker[i];
    }
    // Convert z-fixed boundary attributes into markers
    ess_bdr = 0;
    for (auto zfix_attrib : zfix_attribs)
    {
      ess_bdr[zfix_attrib-1] = 1;
    }
    // Find all z-direction vdofs with the given boundary marker
    par_fe_space.GetEssentialVDofs(ess_bdr, new_ess_vdof_marker, 2);
    // Compute union of x-direction, y-direction, and z-direction vdofs
    for (int i{0}; i < ess_vdof_marker.Size(); ++i)
    {
      ess_vdof_marker[i] = ess_vdof_marker[i] || new_ess_vdof_marker[i];
    }
    // Convert the vdofs to tdofs to remove duplicate values over ranks
    mfem::Array<int> ess_tdof_marker;
    par_fe_space.GetRestrictionMatrix()->BooleanMult(ess_vdof_marker, ess_tdof_marker);
    // Convert the tdof marker array to a tdof list
    mfem::FiniteElementSpace::MarkerToList(ess_tdof_marker, ess_tdof_list);
  }
  timer.stop();
  SLIC_INFO_ROOT(axom::fmt::format(
    "Time to set up boundary conditions: {0:f}ms", timer.elapsedTimeInMilliSec()
  ));

  // This block of code constructs a small-deformation elasticity bilinear form.
  timer.start();
  mfem::ParBilinearForm a(&par_fe_space);
  mfem::ConstantCoefficient lambda_coeff(lambda);
  mfem::ConstantCoefficient mu_coeff(mu);
  a.AddDomainIntegrator(new mfem::ElasticityIntegrator(lambda_coeff, mu_coeff));

  // Assemble the on-rank bilinear form stiffness matrix.
  a.Assemble();
  // Reduce to tdofs and form a hypre parallel matrix for parallel solution of
  // the linear system.
  auto A = std::make_unique<mfem::HypreParMatrix>();
  a.FormSystemMatrix(ess_tdof_list, *A);
  timer.stop();
  SLIC_INFO_ROOT(axom::fmt::format(
    "Time to create and assemble internal stiffness: {0:f}ms", timer.elapsedTimeInMilliSec()
  ));

  // This block of code does initial setup of Tribol.
  timer.start();
  // First, Tribol is initialized with the dimension of the space and the MPI
  // communicator. These are stored globally.
  tribol::initialize(pmesh->SpaceDimension(), MPI_COMM_WORLD);
  // Next, we create a Tribol coupling scheme between the contact surfaces on
  // the MFEM mesh. To create the coupling scheme requires several steps: 1)
  // building a boundary submesh, 2) building a LOR mesh (if required), 3)
  // re-decomposing the domain to move spatially close surface element pairs on
  // to the same rank, and 4) creating a Tribol mesh. Steps 1 and 2 are
  // performed when this method is called. Steps 3 and 4 are accomplished via a
  // call to updateMfemParallelDecomposition(), which is typically called before
  // calling update().
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
  // The basic Lagrange multiplier options are set here. For this problem, we
  // ask Tribol to compute a contact residual and a Jacobian (though we only use
  // the Jacobian).
  tribol::setLagrangeMultiplierOptions(
    coupling_scheme_id,
    tribol::ImplicitEvalMode::MORTAR_RESIDUAL_JACOBIAN
  );

  // save initial configuration
  paraview_datacoll.Save();
  visit_datacoll.Save();
  
  // Update the cycle information for the data collections.
  paraview_datacoll.SetCycle(1);
  paraview_datacoll.SetTime(1.0);
  paraview_datacoll.SetTimeStep(1.0);
  visit_datacoll.SetCycle(1);
  visit_datacoll.SetTime(1.0);
  visit_datacoll.SetTimeStep(1.0);

  // This creates the parallel adjacency-based mesh redecomposition. It also
  // constructs new Tribol meshes as subdomains of the redecomposed mesh.
  tribol::updateMfemParallelDecomposition();
  double dt {1.0};  // time is arbitrary here (no timesteps)
  // This API call computes the contact response and Jacobian given the current
  // mesh configuration.
  tribol::update(1, 1.0, dt);

  // HypreParMatrices holding the Jacobian from contact and stored in an MFEM
  // block operator are returned from this API call.
  auto A_blk = tribol::getMfemBlockJacobian(coupling_scheme_id);
  // Add the Jacobian from the elasticity bilinear form to the top left block
  A_blk->SetBlock(0, 0, A.release());
  timer.stop();
  SLIC_INFO_ROOT(axom::fmt::format(
    "Time to setup Tribol and compute Jacobian: {0:f}ms", timer.elapsedTimeInMilliSec()
  ));

  timer.start();
  // Create a block RHS vector storing forces and gaps. Note no external forces
  // are present in this problem.
  mfem::BlockVector B_blk { A_blk->RowOffsets() };
  B_blk = 0.0;
  // Create a block solution vector storing displacement and pressures.
  mfem::BlockVector X_blk { A_blk->ColOffsets() };
  X_blk = 0.0;

  // This API call saves nodal mortar gap normal contributions to an
  // (uninitialized) grid function. The function sizes and initializes the grid
  // function.
  mfem::Vector g;
  tribol::getMfemGap(coupling_scheme_id, g);

  // Apply a restriction operator on the submesh: maps dofs stored in g to tdofs
  // stored in G for parallel solution of the linear system.
  {
    auto& G = B_blk.GetBlock(1);
    auto& R_submesh = *tribol::getMfemPressure(coupling_scheme_id)
      .ParFESpace()->GetRestrictionOperator();
    R_submesh.Mult(g, G);
  }

  // Use a linear solver to find the block displacement/pressure vector.
  mfem::MINRESSolver solver(MPI_COMM_WORLD);
  solver.SetRelTol(1.0e-8);
  solver.SetAbsTol(1.0e-12);
  solver.SetMaxIter(5000);
  solver.SetPrintLevel(3);
  solver.SetOperator(*A_blk);
  solver.Mult(B_blk, X_blk);

  // Move the block displacements to the displacement grid function.
  {
    auto& U = X_blk.GetBlock(0);
    auto& P = *par_fe_space.GetProlongationMatrix();
    P.Mult(U, displacement);
  }
  // Fix the sign on the displacements.
  displacement.Neg();

  // Update mesh coordinates given the displacement.
  coords += displacement;
  timer.stop();
  SLIC_INFO_ROOT(axom::fmt::format(
    "Time to solve for updated displacements: {0:f}ms", timer.elapsedTimeInMilliSec()
  ));

  // Save the deformed configuration
  paraview_datacoll.Save();
  visit_datacoll.Save();

  // Tribol cleanup: deletes the coupling schemes and clears associated memory.
  tribol::finalize();
  MPI_Finalize();

  return 0;
}
