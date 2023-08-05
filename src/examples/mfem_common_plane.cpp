// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Tribol Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (MIT)

/**
 * @file mfem_common_plane.cpp
 *
 * @brief Demonstrates explicit contact using Tribol's common plane algorithm
 *
 * Demonstrates contact in an explicit finite element code through a simple two
 * block impact problem in three dimensions. The first block occupies [0, 1]^3
 * and the second block occupies [0, 1]x[0, 1]x[1.01, 2.01]. An initial velocity
 * is applied to the second block in the negative z-direction and fixed velocity
 * boundary conditions are applied to the first block at z = 0. Contact
 * constraints are enforced using the common plane algorithm.
 *
 * The example uses the Tribol MFEM interface, which supports decomposed (MPI)
 * meshes and has experimental support for higher order meshes using low-order
 * refinement of higher-order geometry representations. Comments in the main
 * function below give details on each step of the example code.
 *
 * Example runs (from repo root directory):
 *   - mpirun -np 4 {build_dir}/examples/mfem_common_plane_ex
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
#include "tribol/utils/TestUtils.hpp"

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
  // kinematic penalty parameter
  double p_kine = 10000.0;
  // number of cycles to skip before next output
  int output_cycles = 20;
  // material density
  double rho = 1000.0;
  // Lame parameter lambda
  double lambda = 100000.0;
  // Lame parameter mu (shear modulus)
  double mu = 100000.0;

  // parse command line options
  axom::CLI::App app { "mfem_common_plane" };
  app.add_option("-r,--refine", ref_levels,
    "Number of times to refine the mesh uniformly.")
    ->capture_default_str();
  app.add_option("-o,--order", order, 
    "Finite element order (polynomial degree).")
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
  app.add_option("-p,--kinematicpenalty", p_kine, 
    "Kinematic penalty parameter.")
    ->capture_default_str();
  app.add_option("-c,--outputcycles", output_cycles, 
    "Number of cycles to skip before next output.")
    ->capture_default_str();
  app.add_option("-R,--rho", rho, 
    "Material density.")
    ->capture_default_str();
  app.add_option("-l,--lambda", lambda, 
    "Lame parameter lambda.")
    ->capture_default_str();
  app.add_option("-m,--mu", mu, 
    "Lame parameter mu (shear modulus).")
    ->capture_default_str();
  CLI11_PARSE(app, argc, argv);

  SLIC_INFO_ROOT("Running mfem_common_plane with the following options:");
  SLIC_INFO_ROOT(axom::fmt::format("refinement levels: {0}", ref_levels));
  SLIC_INFO_ROOT(axom::fmt::format("polynomial order:  {0}", order));
  SLIC_INFO_ROOT(axom::fmt::format("initial velocity:  {0}", init_velocity));
  SLIC_INFO_ROOT(axom::fmt::format("timestep size:     {0}", dt));
  SLIC_INFO_ROOT(axom::fmt::format("end time:          {0}", t_end));
  SLIC_INFO_ROOT(axom::fmt::format("kinematic penalty: {0}", p_kine));
  SLIC_INFO_ROOT(axom::fmt::format("output interval:   {0}", output_cycles));
  SLIC_INFO_ROOT(axom::fmt::format("density:           {0}", rho));
  SLIC_INFO_ROOT(axom::fmt::format("lambda:            {0}", lambda));
  SLIC_INFO_ROOT(axom::fmt::format("mu:                {0}\n", mu));

  // fixed options
  // location of mesh file. TRIBOL_REPO_DIR is defined in tribol/config.hpp
  std::string mesh_file = TRIBOL_REPO_DIR "/data/two_hex_apart.mesh";
  // boundary element attributes of contact surface 1, the z = 1 plane of the
  // first block
  auto contact_surf_1 = std::set<int>({4});
  // boundary element attributes of contact surface 2, the z = 1.01 plane of the
  // second block
  auto contact_surf_2 = std::set<int>({5});
  // boundary element attributes of fixed surface (z = 0 plane of the first
  // block, all t)
  auto fixed_attrs = std::set<int>({3});
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
  auto paraview_datacoll = mfem::ParaViewDataCollection("common_plane_pv", pmesh.get());
  auto visit_datacoll = mfem::VisItDataCollection("common_plane_vi", pmesh.get());

  // This block of code creates position, displacement, and velocity grid
  // functions (and associated finite element collections and finite element
  // spaces) on the mesh. The displacement and velocity grid functions are
  // initialized to zero. The position grid function is initialized with the
  // nodal coordinates. These grid functions are registered with the data
  // collections for output.
  timer.start();
  // Finite element collection (shared between all grid functions).
  mfem::H1_FECollection fe_coll { order, pmesh->SpaceDimension() };
  // Finite element space (shared between all grid functions).
  mfem::ParFiniteElementSpace par_fe_space {
    pmesh.get(), &fe_coll, pmesh->SpaceDimension() };
  // Create coordinate grid function
  mfem::ParGridFunction coords { &par_fe_space };
  // Set coordinate grid function based on nodal locations. In MFEM, nodal
  // locations of higher order meshes are stored in a grid function. For linear
  // MFEM meshes, nodal locations can be stored in a grid function or through
  // the vertex coordinates. Since we track coordinates in a grid function even
  // for linear meshes, it makes sense to create a nodes grid function so we
  // don't have to manually update the vertices.
  pmesh->SetNodalGridFunction(&coords, false);
  paraview_datacoll.RegisterField("pos", &coords);
  visit_datacoll.RegisterField("pos", &coords);

  // Save reference coordinates
  mfem::ParGridFunction ref_coords { coords };

  // Create a grid function for displacement
  mfem::ParGridFunction displacement { &par_fe_space };
  paraview_datacoll.RegisterField("disp", &displacement);
  visit_datacoll.RegisterField("disp", &displacement);
  displacement = 0.0;

  // Create a grid function for velocity
  mfem::ParGridFunction velocity { &par_fe_space };
  paraview_datacoll.RegisterField("vel", &velocity);
  visit_datacoll.RegisterField("vel", &velocity);
  velocity = 0.0;
  timer.stop();
  SLIC_INFO_ROOT(axom::fmt::format(
    "Time to create grid functions: {0:f}ms", timer.elapsedTimeInMilliSec()
  ));

  timer.start();
  {
    // This block of code sets initial conditions on the velocity grid function.
    // We create a piecewise constant vector coefficient with an initial
    // velocity in the z-component of the second block and zero elsewhere.
    // First, we create a constant vector coefficient, then create arrays of the
    // vector coefficient and the element attributes which belong to the second
    // block to build the piecewise constant vector coefficient. In the
    // PWVectorCoefficient class, zeros are assigned to missing attribute
    // numbers. Finally, the velocity coefficient is projected onto the velocity
    // grid function.
    mfem::Vector init_velocity_vector({0.0, 0.0, -std::abs(init_velocity)});
    mfem::VectorConstantCoefficient init_velocity_coeff(init_velocity_vector);
    mfem::Array<int> moving_attrs_array;
    mfem::Array<mfem::VectorCoefficient*> init_velocity_coeff_array;
    moving_attrs_array.Reserve(moving_attrs.size());
    init_velocity_coeff_array.Reserve(moving_attrs.size());
    for (auto moving_attr : moving_attrs)
    {
      moving_attrs_array.Append(moving_attr);
      init_velocity_coeff_array.Append(&init_velocity_coeff);
    }
    mfem::PWVectorCoefficient initial_v_coeff(
      pmesh->SpaceDimension(), 
      moving_attrs_array,
      init_velocity_coeff_array
    );
    velocity.ProjectCoefficient(initial_v_coeff);
  }

  // This block of code builds a list of degrees of freedom in the x, y, and z
  // directions in the first block along the z = 0 plane. A homogeneous velocity
  // boundary condition is applied to these degrees of freedom.
  mfem::Array<int> ess_vdof_list;
  {
    // First, build an array of "markers" (i.e. booleans) to denote which vdofs
    // are in the list.
    mfem::Array<int> ess_vdof_marker;
    // Also, convert active boundary attributes into markers.
    mfem::Array<int> ess_bdr(pmesh->bdr_attributes.Max());
    ess_bdr = 0;
    for (auto fixed_attr : fixed_attrs)
    {
      ess_bdr[fixed_attr-1] = 1;
    }
    // Note no component is given as an argument here, so all components (x, y,
    // and z) are returned.
    par_fe_space.GetEssentialVDofs(ess_bdr, ess_vdof_marker);
    // Convert the marker array to a list.
    mfem::FiniteElementSpace::MarkerToList(ess_vdof_marker, ess_vdof_list);
  }
  timer.stop();
  SLIC_INFO_ROOT(axom::fmt::format(
    "Time to set up boundary conditions: {0:f}ms", timer.elapsedTimeInMilliSec()
  ));

  // This block of code builds a small-deformation elasticity explicit, lumped
  // mass update operator. Lumping is performed using row summation.
  timer.start();
  mfem::ConstantCoefficient rho_coeff {rho};
  mfem::ConstantCoefficient lambda_coeff {lambda};
  mfem::ConstantCoefficient mu_coeff {mu};
  mfem_ext::ExplicitMechanics op {par_fe_space, rho_coeff, lambda_coeff, mu_coeff};
  timer.stop();
  SLIC_INFO_ROOT(axom::fmt::format(
    "Time to set up elasticity bilinear form: {0:f}ms", timer.elapsedTimeInMilliSec()
  ));

  // This block of code sets up a central difference time integration scheme.
  // The constructor takes a list of degrees of freedom in the velocity grid
  // function that will have homogeneous boundary conditions applied to them.
  mfem_ext::CentralDiffSolver solver { ess_vdof_list };
  solver.Init(op);

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
    *pmesh, coords, contact_surf_1, contact_surf_2,
    tribol::SURFACE_TO_SURFACE, 
    tribol::NO_CASE, 
    tribol::COMMON_PLANE, 
    tribol::FRICTIONLESS,
    tribol::PENALTY,
    tribol::BINNING_GRID
  );
  // This API call adds a velocity field to the coupling scheme. This is used
  // for computing the maximum common plane timestep and, if activated, gap rate
  // penalty.
  tribol::registerMfemVelocity(coupling_scheme_id, velocity);
  // The type of penalty enforcement is set here (i.e. rate vs. kinematic and
  // the method of setting the penalty value).
  tribol::setPenaltyOptions(
    coupling_scheme_id,
    tribol::KINEMATIC,
    tribol::KINEMATIC_CONSTANT,
    tribol::NO_RATE_PENALTY
  );
  timer.stop();
  SLIC_INFO_ROOT(axom::fmt::format(
    "Time to set up Tribol: {0:f}ms", timer.elapsedTimeInMilliSec()
  ));

  // This is the main timestepping loop.
  int cycle {0};
  for (double t {0.0}; t < t_end; t+=dt)
  {
    // Update the cycle information for the data collections.
    paraview_datacoll.SetCycle(cycle);
    paraview_datacoll.SetTime(t);
    paraview_datacoll.SetTimeStep(dt);
    visit_datacoll.SetCycle(cycle);
    visit_datacoll.SetTime(t);
    visit_datacoll.SetTimeStep(dt);

    // Write output if we are at the right cycle.
    if (cycle % output_cycles == 0)
    {
      paraview_datacoll.Save();
      visit_datacoll.Save();
    }

    // This updates the parallel adjacency-based mesh redecomposition based on
    // the latest coordinates. It also constructs new Tribol meshes as
    // subdomains of the redecomposed mesh.
    tribol::updateMfemParallelDecomposition();
    // The next two API calls set the penalty data on the newly constructed
    // Tribol meshes.
    tribol::setKinematicConstantPenalty(mesh1_id, p_kine);
    tribol::setKinematicConstantPenalty(mesh2_id, p_kine);
    // This API call computes the contact response given the current mesh
    // configuration.
    tribol::update(cycle, t, dt);
    // Nodal forces are returned to an MFEM grid function using this API
    // function. We set the contact nodal forces equal to the external forces in
    // the explicit mechanics update operator.
    op.f_ext = 0.0;
    tribol::getMfemResponse(coupling_scheme_id, op.f_ext);

    // Update the time in the explicit mechanics update.
    op.SetTime(t);
    // Take a step in time in the time integration scheme.
    solver.Step(displacement, velocity, t, dt);

    // Update the coordinates based on the new displacement.
    coords.Set(1.0, ref_coords);
    coords += displacement;

    // Increment the cycle count.
    ++cycle;
  }

  // save output after last timestep
  paraview_datacoll.Save();
  visit_datacoll.Save();

  // Tribol cleanup: deletes the coupling schemes and clears associated memory.
  tribol::finalize();
  MPI_Finalize();

  return 0;
}
