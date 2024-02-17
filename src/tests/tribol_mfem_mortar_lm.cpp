// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Tribol Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (MIT)

#include <set>

#include <gtest/gtest.h>

// Tribol includes
#include "tribol/config.hpp"
#include "tribol/common/Parameters.hpp"
#include "tribol/interface/tribol.hpp"
#include "tribol/interface/mfem_tribol.hpp"
#include "tribol/types.hpp"

// Redecomp includes
#include "redecomp/redecomp.hpp"

#ifdef TRIBOL_USE_UMPIRE
// Umpire includes
#include "umpire/ResourceManager.hpp"
#endif

// MFEM includes
#include "mfem.hpp"

// Axom includes
#include "axom/CLI11.hpp"
#include "axom/slic.hpp"

/**
 * @brief This tests the Tribol MFEM interface running a contact patch test.
 *
 */
class MfemCommonPlaneTest : public testing::TestWithParam<int> {
protected:
  double max_disp_;
  void SetUp() override
  {
    // number of times to uniformly refine the serial mesh before constructing the
    // parallel mesh
    int ref_levels = GetParam();
    // polynomial order of the finite element discretization
    int order = 1;

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
    // boundary element attributes of z-fixed surfaces (3: surface at z = 0, 6: surface at z = 1.99)
    auto zfixed_attrs = std::set<int>({3, 6});

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

    // grid function for displacement
    mfem::ParGridFunction displacement { &par_fe_space };
    displacement = 0.0;

    // recover dirichlet bc tdof list
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
      for (auto zfixed_attr : zfixed_attrs)
      {
        ess_bdr[zfixed_attr-1] = 1;
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
    mfem::ParBilinearForm a(&par_fe_space);
    mfem::ConstantCoefficient lambda(50.0);
    mfem::ConstantCoefficient mu(50.0);
    a.AddDomainIntegrator(new mfem::ElasticityIntegrator(lambda, mu));
    a.Assemble();

    // compute elasticity contribution to stiffness
    auto A = std::make_unique<mfem::HypreParMatrix>();
    a.FormSystemMatrix(ess_tdof_list, *A);

    // set up tribol
    tribol::initialize(pmesh->SpaceDimension(), MPI_COMM_WORLD);
    coords.ReadWrite();
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

    coords.ReadWrite();
    // update tribol (compute contact contribution to force and stiffness)
    tribol::updateMfemParallelDecomposition();
    double dt {1.0};  // time is arbitrary here (no timesteps)
    tribol::update(1, 1.0, dt);

    // retrieve block stiffness matrix
    auto A_blk = tribol::getMfemBlockJacobian(0);
    A_blk->SetBlock(0, 0, A.release());

    // create block solution and RHS vectors
    mfem::Vector B_blk { A_blk->Height() };
    B_blk.UseDevice(true);
    B_blk = 0.0;
    mfem::Vector X_blk { A_blk->Width() };
    X_blk.UseDevice(true);
    X_blk = 0.0;

    // retrieve gap vector (RHS) from contact
    mfem::Vector g;
    g.UseDevice(true);
    tribol::getMfemGap(0, g);

    // prolongation transpose operator on submesh: maps dofs stored in g to tdofs stored in G
    int n_disp_dofs = par_fe_space.GetTrueVSize();
    auto& pressure = tribol::getMfemPressure(coupling_scheme_id);
    int n_lm_dofs = pressure.ParFESpace()->GetTrueVSize();
    {
      mfem::Vector G(B_blk, n_disp_dofs, n_lm_dofs);
      auto& P_submesh = *pressure.ParFESpace()->GetProlongationMatrix();
      P_submesh.MultTranspose(g, G);
      B_blk.SyncMemory(G);
    }

    mfem::Array2D<mfem::HypreParMatrix*> hypre_blocks(2, 2);
    for (int i{0}; i < 2; ++i)
    {
      for (int j{0}; j < 2; ++j)
      {
        if (A_blk->GetBlock(i, j).Height() != 0 && A_blk->GetBlock(i, j).Width() != 0)
        {
          hypre_blocks(i, j) = const_cast<mfem::HypreParMatrix*>(
            dynamic_cast<const mfem::HypreParMatrix*>(&A_blk->GetBlock(i, j))
          );
        }
        else
        {
          hypre_blocks(i, j) = nullptr;
        }
      }
    }
    auto A_merged = std::unique_ptr<mfem::HypreParMatrix>(
      mfem::HypreParMatrixFromBlocks(hypre_blocks)
    );

    // solve for X_blk
    mfem::MINRESSolver solver(MPI_COMM_WORLD);
    solver.SetRelTol(1.0e-8);
    solver.SetAbsTol(1.0e-12);
    solver.SetMaxIter(5000);
    solver.SetPrintLevel(3);
    solver.SetOperator(*A_merged);
    mfem::HypreDiagScale prec(*A_merged);
    solver.SetPreconditioner(prec);
    solver.Mult(B_blk, X_blk);

    // move block displacements to grid function
    {
      mfem::Vector U(X_blk, 0, n_disp_dofs);
      auto& P = *par_fe_space.GetProlongationMatrix();
      P.Mult(U, displacement);
    }
    displacement.Neg();

    auto local_max = displacement.Max();
    max_disp_ = 0.0;
    MPI_Allreduce(&local_max, &max_disp_, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  }
};

TEST_P(MfemCommonPlaneTest, mass_matrix_transfer)
{
  EXPECT_LT(std::abs(max_disp_ - 0.005), 1.0e-6);

  MPI_Barrier(MPI_COMM_WORLD);
}

INSTANTIATE_TEST_SUITE_P(tribol, MfemCommonPlaneTest, testing::Values(2));

//------------------------------------------------------------------------------
#include "axom/slic/core/SimpleLogger.hpp"

int main(int argc, char* argv[])
{
  int result = 0;

  MPI_Init(&argc, &argv);

  ::testing::InitGoogleTest(&argc, argv);

#ifdef TRIBOL_USE_UMPIRE
  umpire::ResourceManager::getInstance();  // initialize umpire's ResouceManager
#endif

  axom::slic::SimpleLogger logger;  // create & initialize test logger, finalized when
                                    // exiting main scope

#ifdef TRIBOL_ENABLE_CUDA
  mfem::Device device("cuda");
#endif

  result = RUN_ALL_TESTS();

  tribol::finalize();
  MPI_Finalize();

  return result;
}
