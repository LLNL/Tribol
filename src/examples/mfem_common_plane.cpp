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
 * block impact problem.  Contact constraints are enforced using the common
 * plane algorithm.
 *
 * The example uses the Tribol MFEM interface, which supports decomposed (MPI)
 * meshes and has experimental support for higher order meshes using low-order
 * refinement of higher-order geometry representations.
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

  // command line options
  // number of times to uniformly refine the serial mesh before constructing the
  // parallel mesh
  int ref_levels = 2;
  // polynomial order of the finite element discretization
  int order = 1;
  // initial velocity
  double initial_v = 0.001;
  // timestep size
  double dt = 0.001;
  // end time
  double t_end = 0.35;
  // kinematic penalty
  double p_kine = 500.0;
  // number of cycles to skip before output
  int output_cycles = 5;
  // material density
  double rho = 100.0;
  // lame parameter
  double lambda = 100000.0;
  // lame parameter (shear modulus)
  double mu = 100000.0;

  axom::CLI::App app { "mfem_common_plane" };
  app.add_option("-r,--refine", ref_levels,
    "Number of times to refine the mesh uniformly.")
    ->capture_default_str();
  app.add_option("-o,--order", order, 
    "Finite element order (polynomial degree).")
    ->capture_default_str();
  app.add_option("-v,--initialv", initial_v, 
    "Initial velocity of the top block.")
    ->capture_default_str();
  app.add_option("-d,--dt", dt, 
    "Timestep size.")
    ->capture_default_str();
  app.add_option("-e,--endtime", t_end, 
    "Time of the end of the simulation.")
    ->capture_default_str();
  app.add_option("-p,--kinematicpenalty", p_kine, 
    "Kinematic penalty parameter.")
    ->capture_default_str();
  app.add_option("-c,--outputcycles", output_cycles, 
    "Cycles to skip before next output.")
    ->capture_default_str();
  app.add_option("-R,--rho", rho, 
    "Material density.")
    ->capture_default_str();
  app.add_option("-l,--lambda", lambda, 
    "Lame parameter.")
    ->capture_default_str();
  app.add_option("-m,--mu", mu, 
    "Lame parameter (shear modulus).")
    ->capture_default_str();
  CLI11_PARSE(app, argc, argv);

  SLIC_INFO_ROOT("Running mfem_common_plane with the following options:");
  SLIC_INFO_ROOT(axom::fmt::format("refine: {0}", ref_levels));
  SLIC_INFO_ROOT(axom::fmt::format("order:  {0}\n", order));

  // fixed options
  // location of mesh file. TRIBOL_REPO_DIR is defined in tribol/config.hpp
  std::string mesh_file = TRIBOL_REPO_DIR "/data/two_hex_apart.mesh";
  // boundary element attributes of contact surface 1
  auto contact_surf_1 = std::set<int>({4});
  // boundary element attributes of contact surface 2
  auto contact_surf_2 = std::set<int>({5});
  // boundary element attributes of fixed surface (points on z = 0, all t)
  auto fixed_attrs = std::set<int>({3});
  // element attribute corresponding to volume elements where an initial
  // velocity will be applied
  auto moving_attrs = std::set<int>({2});

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
  auto paraview_datacoll = mfem::ParaViewDataCollection("common_plane_pv", pmesh.get());
  auto visit_datacoll = mfem::VisItDataCollection("common_plane_vi", pmesh.get());

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

  // grid function for velocity
  mfem::ParGridFunction velocity { &par_fe_space };
  paraview_datacoll.RegisterField("vel", &velocity);
  visit_datacoll.RegisterField("vel", &velocity);
  velocity = 0.0;
  timer.stop();
  SLIC_INFO_ROOT(axom::fmt::format(
    "Time to create grid functions: {0:f}ms", timer.elapsedTimeInMilliSec()
  ));

  // set initial velocity
  timer.start();
  {
    mfem::Array<int> attrib_marker(pmesh->attributes.Max());
    attrib_marker = 0;
    for (auto moving_attr : moving_attrs)
    {
      attrib_marker[moving_attr-1] = 1;
    }
    mfem::Array<int> vdof_marker(par_fe_space.GetVSize());
    vdof_marker = 0;
    for (int e{0}; e < pmesh->GetNE(); ++e)
    {
      if (attrib_marker[pmesh->GetElement(e)->GetAttribute()-1])
      {
        mfem::Array<int> vdofs;
        par_fe_space.GetElementDofs(e, vdofs);
        par_fe_space.DofsToVDofs(2, vdofs);
        for (auto vdof : vdofs)
        {
          vdof_marker[vdof] = 1;
        }
      }
    }
    for (int i{0}; i < vdof_marker.Size(); ++i)
    {
      if (vdof_marker[i])
      {
        velocity[i] = -std::abs(initial_v);
      }
    }
  }

  // recover dirichlet bc tdof list
  mfem::Array<int> ess_vdof_list;
  {
    mfem::Array<int> ess_vdof_marker;
    mfem::Array<int> ess_bdr(pmesh->bdr_attributes.Max());
    ess_bdr = 0;
    for (auto fixed_attr : fixed_attrs)
    {
      ess_bdr[fixed_attr-1] = 1;
    }
    par_fe_space.GetEssentialVDofs(ess_bdr, ess_vdof_marker);
    mfem::FiniteElementSpace::MarkerToList(ess_vdof_marker, ess_vdof_list);
  }
  timer.stop();
  SLIC_INFO_ROOT(axom::fmt::format(
    "Time to set up boundary conditions: {0:f}ms", timer.elapsedTimeInMilliSec()
  ));

  // set up mfem elasticity bilinear form
  timer.start();
  mfem::ConstantCoefficient rho_coeff {rho};
  mfem::ConstantCoefficient lambda_coeff {lambda};
  mfem::ConstantCoefficient mu_coeff {mu};
  mfem_ext::ExplicitMechanics op {par_fe_space, rho_coeff, lambda_coeff, mu_coeff};
  timer.stop();
  SLIC_INFO_ROOT(axom::fmt::format(
    "Time to set up elasticity bilinear form: {0:f}ms", timer.elapsedTimeInMilliSec()
  ));

  // set up time integrator
  mfem_ext::CentralDiffSolver solver { ess_vdof_list };
  solver.Init(op);

  // set up tribol
  timer.start();
  tribol::initialize(pmesh->SpaceDimension(), MPI_COMM_WORLD);
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
  tribol::registerMfemVelocity(0, velocity);
  tribol::setPenaltyOptions(
    0,
    tribol::KINEMATIC,
    tribol::KINEMATIC_CONSTANT,
    tribol::NO_RATE_PENALTY
  );
  timer.stop();
  SLIC_INFO_ROOT(axom::fmt::format(
    "Time to set up Tribol: {0:f}ms", timer.elapsedTimeInMilliSec()
  ));

  int cycle {0};
  for (double t {0.0}; t < t_end; t+=dt)
  {
    paraview_datacoll.SetCycle(cycle);
    paraview_datacoll.SetTime(t);
    paraview_datacoll.SetTimeStep(dt);
    visit_datacoll.SetCycle(cycle);
    visit_datacoll.SetTime(t);
    visit_datacoll.SetTimeStep(dt);

    if (cycle % output_cycles == 0)
    {
      paraview_datacoll.Save();
      visit_datacoll.Save();
    }

    tribol::updateMfemParallelDecomposition();
    tribol::setKinematicConstantPenalty(0, p_kine);
    tribol::setKinematicConstantPenalty(1, p_kine);
    tribol::update(cycle, t, dt);
    op.f_ext = 0.0;
    tribol::getMfemResponse(0, op.f_ext);

    op.SetTime(t);
    solver.Step(displacement, velocity, t, dt);

    coords += displacement;
    // nodal coordinates stored in Nodes grid function for higher order
    if (order == 1)
    {
      pmesh->SetVertices(coords);
    }

    ++cycle;
  }

  // save output after last timestep
  paraview_datacoll.Save();
  visit_datacoll.Save();

  // cleanup
  tribol::finalize();
  MPI_Finalize();

  return 0;
}
