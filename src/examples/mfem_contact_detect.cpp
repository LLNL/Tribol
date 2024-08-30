// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Tribol Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (MIT)

/**
 * @file mfem_contact_detect.cpp
 *
 * @brief Demonstrates contact detection using the MFEM interface
 *
 * Demonstrates a three dimensional contact detection using the MFEM interface
 *
 * Example runs (from repo root directory):
 *   - mpirun -np 4 {build_dir}/examples/mfem_contact_detect_ex
 *
 * Example output can be viewed in VisIt or ParaView.
 */

#include <set>

// Tribol includes
#include "tribol/config.hpp"
#include "tribol/interface/tribol.hpp"
#include "tribol/interface/mfem_tribol.hpp"

#ifdef TRIBOL_USE_UMPIRE
// Umpire includes
#include "umpire/ResourceManager.hpp"
#endif

// Axom includes
#include "axom/CLI11.hpp"
#include "axom/slic.hpp"

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
  // initial velocity to apply in the negative z-direction to the top block
  double init_velocity = 0.02;
  // timestep size (fixed)
  double dt = 0.001;
  // end time of the simulation
  double t_end = 2.0;
  // number of cycles to skip before next output
  int output_cycles = 20;

  // parse command line options
  axom::CLI::App app { "mfem_contact_detect" };
  app.add_option("-r,--refine", ref_levels,
    "Number of times to refine the mesh uniformly.")
    ->capture_default_str();
  app.add_option("-v,--initialv", init_velocity, 
    "Initial velocity of the top block.")
    ->capture_default_str();
  app.add_option("-d,--dt", dt, 
    "Timestep size (fixed).")
    ->capture_default_str();
  app.add_option("-e,--endtime", t_end, 
    "End time of the simulation.")
    ->capture_default_str();
  app.add_option("-c,--outputcycles", output_cycles, 
    "Number of cycles to skip before next output.")
    ->capture_default_str();
  // TODO: LOR support for implicit contact
  // app.add_option("-o,--order", order, 
  //   "Finite element order (polynomial degree).")
  //   ->capture_default_str();
  CLI11_PARSE(app, argc, argv);

  SLIC_INFO_ROOT("Running mfem_contact_detect with the following options:");
  SLIC_INFO_ROOT(axom::fmt::format("refine: {0}", ref_levels));
  SLIC_INFO_ROOT(axom::fmt::format("initial velocity:  {0}", init_velocity));
  SLIC_INFO_ROOT(axom::fmt::format("timestep size:     {0}", dt));
  SLIC_INFO_ROOT(axom::fmt::format("end time:          {0}", t_end));

  // fixed options
  // location of mesh file. TRIBOL_REPO_DIR is defined in tribol/config.hpp
  std::string mesh_file = TRIBOL_REPO_DIR "/data/two_hex_apart.mesh";
  // boundary element attributes of mortar surface, the z = 1 plane of the first
  // block
  auto mortar_attrs = std::set<int>({4});
  // boundary element attributes of nonmortar surface, the z = 1.01 plane of the
  // second block
  auto nonmortar_attrs = std::set<int>({5});
  // element attribute corresponding to volume elements where an initial
  // velocity will be applied in the z-direction
  auto moving_attrs = std::set<int>({2});
  
  // create an axom timer to give wall times for each step
  axom::utilities::Timer timer { false };

  // This block of code will read the mesh data given in two_hex_apart.mesh,
  // create an mfem::Mesh, refine the mesh, then create an mfem::ParMesh.
  // Optionally, the mfem::ParMesh can be refined further on each rank by
  // setting par_ref_levels >= 1, though this is disabled below.
  timer.start();
  std::unique_ptr<mfem::ParMesh> pmesh { nullptr };
  {
    // read serial mesh
    auto mesh = std::make_unique<mfem::Mesh>(mesh_file.c_str(), 1, 1);

    // refine serial mesh
    for (int i{0}; i < ref_levels; ++i)
    {
      mesh->UniformRefinement();
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

  // This block of code creates position and delta position (i.e. change in
  // position for each timestep) grid functions (and associated finite element
  // collections and finite element spaces) on the mesh. The displacement grid
  // function is initialized to zero. The position grid function is initialized
  // with the nodal coordinates. These grid functions are registered with the
  // data collections for output.
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

  // Create a grid function for delta coords
  mfem::ParGridFunction d_coords { &par_fe_space };
  d_coords = 0.0;

  mfem::Vector init_d_coords_vector({0.0, 0.0, -std::abs(init_velocity * dt)});
  mfem::VectorConstantCoefficient init_d_coords_coeff(init_d_coords_vector);
  mfem::Array<int> moving_attrs_array;
  mfem::Array<mfem::VectorCoefficient*> init_d_coords_coeff_array;
  moving_attrs_array.Reserve(moving_attrs.size());
  init_d_coords_coeff_array.Reserve(moving_attrs.size());
  for (auto moving_attr : moving_attrs)
  {
    moving_attrs_array.Append(moving_attr);
    init_d_coords_coeff_array.Append(&init_d_coords_coeff);
  }
  mfem::PWVectorCoefficient initial_dc_coeff(
    pmesh->SpaceDimension(), 
    moving_attrs_array,
    init_d_coords_coeff_array
  );
  d_coords.ProjectCoefficient(initial_dc_coeff);
  timer.stop();
  SLIC_INFO_ROOT(axom::fmt::format(
    "Time to create grid functions: {0:f}ms", timer.elapsedTimeInMilliSec()
  ));

  // This block of code does initial setup of Tribol.
  timer.start();
  // First, Tribol is initialized with the spatial dimension and the MPI
  // communicator. These are stored globally.
  tribol::initialize(pmesh->SpaceDimension(), MPI_COMM_WORLD);
  // Next, we create a Tribol coupling scheme between the contact surfaces on
  // the MFEM mesh. To create the coupling scheme requires several steps: 1)
  // building a boundary submesh, 2) building a LOR mesh (if required), 3)
  // re-decomposing the domain to move spatially close surface element pairs on
  // to the same rank, 4) creating Tribol meshes of each surface, and 5)
  // registering the meshes and coupling scheme with Tribol. These 5 steps are
  // performed by calling two methods: 1) registerMfemCouplingScheme() (steps 1
  // and 2) and 2) updateMfemParallelDecomposition() (steps 3, 4, and 5).
  // registerMfemCouplingScheme() is called here and
  // updateMfemParallelDecomposition() is typically called before calling
  // update().
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
  // ask Tribol to compute a mortar gap.
  tribol::setLagrangeMultiplierOptions(
    coupling_scheme_id,
    tribol::ImplicitEvalMode::MORTAR_GAP
  );

  // Setup VisIt data collections
  mfem::ParGridFunction gap_gridfn(tribol::getMfemPressure(coupling_scheme_id).ParFESpace());
  mfem::VisItDataCollection visit_vol_dc("contact-detect-volume", pmesh.get());
  mfem::VisItDataCollection visit_surf_dc("contact-detect-surface",
                                           gap_gridfn.ParFESpace()->GetMesh());
  visit_vol_dc.RegisterField("pos", &coords);
  visit_vol_dc.RegisterField("dpos", &d_coords);
  visit_surf_dc.RegisterField("gap", &gap_gridfn);

  // This is the main timestepping loop.
  int cycle {0};
  double init_dt {dt};
  for (double t {0.0}; t < t_end; t+=dt)
  {
    // Update the cycle information for the data collections.
    visit_vol_dc.SetCycle(cycle);
    visit_vol_dc.SetTime(t);
    visit_vol_dc.SetTimeStep(dt);
    visit_surf_dc.SetCycle(cycle);
    visit_surf_dc.SetTime(t);
    visit_surf_dc.SetTimeStep(dt);

    // Write output if we are at the right cycle.
    if (cycle % output_cycles == 0)
    {
      visit_vol_dc.Save();
      visit_surf_dc.Save();
    }

    // This updates the parallel adjacency-based mesh redecomposition based on
    // the latest coordinates. It also constructs new Tribol meshes as
    // subdomains of the redecomposed mesh.
    tribol::updateMfemParallelDecomposition();
    // This API call computes the contact response given the current mesh
    // configuration.
    tribol::update(cycle, t, dt);
    // Reset dt (in case Tribol changed it)
    dt = init_dt;

    mfem::ParLinearForm gap(gap_gridfn.ParFESpace());
    tribol::getMfemGap(coupling_scheme_id, gap);
    std::unique_ptr<mfem::HypreParVector> gap_true(gap.ParallelAssemble());
    gap_gridfn.SetFromTrueDofs(*gap_true);

    double on_rank_min_gap = gap_gridfn.Min();
    double global_min_gap = 0.0;
    MPI_Reduce(&on_rank_min_gap, &global_min_gap, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    if (global_min_gap < 0.0)
    {
      if (rank == 0)
      {
        std::cout << "Contact detected at t = " << t << std::endl;
      }
      break;
    }
    
    // Update the coordinates.
    coords += d_coords;

    // Increment the cycle count.
    ++cycle;
  }

  // save output after last timestep
  visit_vol_dc.Save();
  visit_surf_dc.Save();

  // Tribol cleanup: deletes the coupling schemes and clears associated memory.
  tribol::finalize();
  MPI_Finalize();

  return 0;
}
