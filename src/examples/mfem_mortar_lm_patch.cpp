// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Tribol Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (MIT)

#include <set>

#include "axom/CLI11.hpp"
#include "axom/slic.hpp"

#include "mfem.hpp"

#include "redecomp/redecomp.hpp"
#include "tribol/common/Parameters.hpp"
#include "tribol/config.hpp"
#include "tribol/interface/tribol.hpp"

template <int NDIMS>
void RedecompExample(mfem::ParMesh& pmesh, int order, double max_out_of_balance);

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
  // number of times to uniformly refine the serial mesh before constructing the
  // parallel mesh
  int ref_levels = 0;
  // polynomial order of the finite element discretization
  int order = 1;

  axom::CLI::App app { "mfem_mortar_lm_patch" };
  app.add_option("-r,--refine", ref_levels,
    "Number of times to refine the mesh uniformly.")
    ->capture_default_str();
  app.add_option("-o,--order", order, 
    "Finite element order (polynomial degree).")
    ->capture_default_str();
  CLI11_PARSE(app, argc, argv);

  SLIC_INFO_ROOT("Running mfem_mortar_lm_patch with the following options:");
  SLIC_INFO_ROOT(axom::fmt::format("refine: {0}", ref_levels));
  SLIC_INFO_ROOT(axom::fmt::format("order:  {0}\n", order));

  // fixed options
  // location of mesh file. TRIBOL_REPO_DIR is defined in tribol/config.hpp
  std::string mesh_file = TRIBOL_REPO_DIR "/data/two_hex_overlap.mesh";
  // boundary element attributes of mortar surface
  auto mortar_attribs = std::set<int>({4});
  // boundary element attributes of nonmortar surface
  auto nonmortar_attribs = std::set<int>({5});
  // boundary element attributes of x-fixed surfaces
  auto xfix_attribs = std::set<int>({1});
  // boundary element attributes of y-fixed surfaces
  auto yfix_attribs = std::set<int>({2});
  // boundary element attributes of z-fixed surfaces
  auto zfix_attribs = std::set<int>({3, 6});

  // read mesh
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
    axom::utilities::Timer timer { false };
    timer.start();
    pmesh = std::make_unique<mfem::ParMesh>(MPI_COMM_WORLD, *mesh);
    mesh.reset(nullptr);
    timer.stop();
    SLIC_INFO_ROOT(axom::fmt::format(
      "Time to create parallel mesh: {0:f}ms", timer.elapsedTimeInMilliSec()
    ));

    // further refinement of parallel mesh
    {
      int par_ref_levels = 2;
      for (int i{0}; i < par_ref_levels; ++i)
      {
        pmesh->UniformRefinement();
      }
    }
  }

  // grid function for higher-order nodes
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

  // recover dirichlet bc tdof list
  mfem::Array<int> ess_tdof_list;
  {
    mfem::Array<int> ess_vdof_marker;
    mfem::Array<int> ess_bdr(pmesh->bdr_attributes.Max());
    ess_bdr = 0;
    for (auto xfix_attrib : xfix_attribs)
    {
      ess_bdr[xfix_attrib] = 1;
    }
    par_fe_space.GetEssentialVDofs(ess_bdr, ess_vdof_marker, 0);
    mfem::Array<int> new_ess_vdof_marker;
    ess_bdr = 0;
    for (auto yfix_attrib : yfix_attribs)
    {
      ess_bdr[yfix_attrib] = 1;
    }
    par_fe_space.GetEssentialVDofs(ess_bdr, new_ess_vdof_marker, 1);
    for (int i{0}; i < ess_vdof_marker.Size(); ++i)
    {
      ess_vdof_marker[i] = ess_vdof_marker[i] || new_ess_vdof_marker[i];
    }
    ess_bdr = 0;
    for (auto zfix_attrib : zfix_attribs)
    {
      ess_bdr[zfix_attrib] = 1;
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

  // set up mfem elasticity bilinear form
  mfem::ConstantCoefficient lambda(50.0);
  mfem::ConstantCoefficient mu(50.0);
  mfem::ParBilinearForm a(&par_fe_space);
  a.AddDomainIntegrator(new mfem::ElasticityIntegrator(lambda, mu));
  a.Assemble();

  // compute elasticity contribution to stiffness
  mfem::HypreParMatrix A;
  a.FormSystemMatrix(ess_tdof_list, A);

  // set up tribol
  tribol::initialize(pmesh->SpaceDimension(), MPI_COMM_WORLD);
  tribol::registerParMesh(
    0, 0, 1, *pmesh, coords, mortar_attribs, nonmortar_attribs,
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
  double dt {1.0};
  tribol::update(1, 1.0, dt);

  // retrieve contact contribution to force
  auto b = tribol::getResponseGridFn(0);

  // retrieve nodal mortar gaps
  auto g = tribol::getGapGridFn(0);

  // retrieve nodal mortar pressures
  auto& p = tribol::getPressureGridFn(0);

  // retrieve block stiffness matrix
  auto A_blk = tribol::getMfemBlockJacobian(0);

  // cleanup
  MPI_Finalize();

  return 0;
}
