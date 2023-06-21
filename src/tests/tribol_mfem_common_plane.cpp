// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Tribol Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (MIT)

#include <set>

#include <gtest/gtest.h>

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


/// Simple central difference method.
class CentralDiffSolver : public mfem::SecondOrderODESolver
{
public:
  CentralDiffSolver(const mfem::Array<int>& bc_vdofs_)
  : bc_vdofs { bc_vdofs_ },
    first_step { true } {}

  void Step(mfem::Vector& x, mfem::Vector& dxdt, double& t, double& dt) override
  {
    // acceleration at t
    f->SetTime(t);
    if (first_step)
    {
      accel.SetSize(x.Size());
      f->Mult(x, dxdt, accel);
      first_step = false;
    }

    // velocity at t + dt/2
    dxdt.Add(0.5*dt, accel);

    // set homogeneous velocity BC at t + dt/2
    SetHomogeneousBC(dxdt);

    // set displacement at t + dt
    x.Add(dt, dxdt);

    // acceleration at t + dt
    f->SetTime(t + dt);
    f->Mult(x, dxdt, accel);

    // velocity at t + dt
    dxdt.Add(0.5*dt, accel);
  }
private:
  mfem::Vector accel;
  mfem::Array<int> bc_vdofs;
  bool first_step;

  void SetHomogeneousBC(mfem::Vector& dxdt) const
  {
    for (auto bc_vdof : bc_vdofs)
    {
      dxdt[bc_vdof] = 0.0;
    }
  }
};

/// Explicit solid mechanics update (lumped mass, optional damping)
class ExplicitMechanics : public mfem::SecondOrderTimeDependentOperator
{
public:
  ExplicitMechanics(
    mfem::ParFiniteElementSpace& fespace, 
    mfem::Coefficient& rho,
    mfem::Coefficient& lambda,
    mfem::Coefficient& mu,
    std::unique_ptr<mfem::ParBilinearForm> damping_ = nullptr
  )
  : elasticity { &fespace },
    damping { std::move(damping_) },
    inv_lumped_mass { fespace.GetVSize() }
  {
    mfem::ParBilinearForm mass { &fespace };
    mass.AddDomainIntegrator(new mfem::VectorMassIntegrator(rho));
    mass.Assemble();
    mfem::Vector ones {fespace.GetVSize()};
    ones = 1.0;
    mass.SpMat().Mult(ones, inv_lumped_mass);
    mfem::Vector mass_true(fespace.GetTrueVSize());
    const Operator& P = *fespace.GetProlongationMatrix();
    P.MultTranspose(inv_lumped_mass, mass_true);
    for (int i {0}; i < mass_true.Size(); ++i)
    {
      mass_true[i] = 1.0 / mass_true[i];
    }
    P.Mult(mass_true, inv_lumped_mass);

    elasticity.AddDomainIntegrator(new mfem::ElasticityIntegrator(lambda, mu));
    elasticity.Assemble();
  }

  void Mult(
    const mfem::Vector& u,
    const mfem::Vector& dudt,
    mfem::Vector& a
  ) const override
  {
    mfem::Vector f { u.Size() };
    f = 0.0;

    mfem::Vector f_int { f.Size() };
    elasticity.Mult(u, f_int);
    f.Add(-1.0, f_int);

    if (damping)
    {
      damping->Mult(dudt, f_int);
      f.Add(-1.0, f_int);
    }

    // sum forces over ranks
    auto& fespace = *elasticity.ParFESpace();
    const Operator& P = *fespace.GetProlongationMatrix();
    mfem::Vector f_true {fespace.GetTrueVSize()};
    P.MultTranspose(f, f_true);
    P.Mult(f_true, f);

    // force already summed over ranks
    f.Add(1.0, tribol::getMfemResponse(0));

    for (int i {0}; i < inv_lumped_mass.Size(); ++i)
    {
      a[i] = inv_lumped_mass[i] * f[i];
    }
  }

private:
  mfem::ParBilinearForm elasticity;
  std::unique_ptr<mfem::ParBilinearForm> damping;
  mfem::Vector inv_lumped_mass;
};

/**
 * @brief This tests the Tribol MFEM interface running a small common plane
 * explicit contact example.
 *
 */
class MfemCommonPlaneTest : public testing::TestWithParam<int> {
protected:
  double max_disp_;
  void SetUp() override
  {
    // number of times to uniformly refine the serial mesh before constructing the
    // parallel mesh
    int ref_levels = 2;
    // polynomial order of the finite element discretization
    int order = GetParam();
    // initial velocity
    double initial_v = 0.001;
    // timestep size
    double dt = 0.001;
    // end time
    double t_end = 0.35;
    // kinematic penalty
    double p_kine = 500.0;

    // fixed options
    // location of mesh file. TRIBOL_REPO_DIR is defined in tribol/config.hpp
    std::string mesh_file = TRIBOL_REPO_DIR "/data/two_hex_apart.mesh";
    // boundary element attributes of surface 1
    auto surf1_attribs = std::set<int>({4});
    // boundary element attributes of surface 2
    auto surf2_attribs = std::set<int>({5});
    // boundary element attributes of fixed surface (all t)
    auto fix_attribs = std::set<int>({3});
    // element attributes of moving volumes (initial condition)
    auto move_attribs = std::set<int>({2});

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
    mfem::ParGridFunction u { &par_fe_space };
    u = 0.0;

    // grid function for velocity
    mfem::ParGridFunction v { &par_fe_space };
    v = 0.0;

    // set initial velocity
    {
      mfem::Array<int> attrib_marker(pmesh->attributes.Max());
      attrib_marker = 0;
      for (auto move_attrib : move_attribs)
      {
        attrib_marker[move_attrib-1] = 1;
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
          v[i] = -std::abs(initial_v);
        }
      }
    }

    // recover dirichlet bc tdof list
    mfem::Array<int> ess_vdof_list;
    {
      mfem::Array<int> ess_vdof_marker;
      mfem::Array<int> ess_bdr(pmesh->bdr_attributes.Max());
      ess_bdr = 0;
      for (auto fix_attrib : fix_attribs)
      {
        ess_bdr[fix_attrib-1] = 1;
      }
      par_fe_space.GetEssentialVDofs(ess_bdr, ess_vdof_marker);
      mfem::FiniteElementSpace::MarkerToList(ess_vdof_marker, ess_vdof_list);
    }

    // set up mfem elasticity bilinear form
    mfem::ConstantCoefficient rho {100.0};
    mfem::ConstantCoefficient lambda {100000.0};
    mfem::ConstantCoefficient mu {100000.0};
    ExplicitMechanics op {par_fe_space, rho, lambda, mu};

    // set up time integrator
    CentralDiffSolver solver { ess_vdof_list };
    solver.Init(op);

    // set up tribol
    tribol::initialize(pmesh->SpaceDimension(), MPI_COMM_WORLD);
    tribol::registerMfemMesh(
      0, 0, 1, *pmesh, coords, surf1_attribs, surf2_attribs,
      tribol::SURFACE_TO_SURFACE, 
      tribol::NO_CASE, 
      tribol::COMMON_PLANE, 
      tribol::FRICTIONLESS,
      tribol::PENALTY,
      tribol::BINNING_GRID
    );
    tribol::registerMfemVelocity(0, v);
    tribol::setPenaltyOptions(
      0,
      tribol::KINEMATIC,
      tribol::KINEMATIC_CONSTANT,
      tribol::NO_RATE_PENALTY
    );
    tribol::setKinematicConstantPenalty(0, p_kine);
    tribol::setKinematicConstantPenalty(1, p_kine);

    int cycle {0};
    for (double t {0.0}; t < t_end; t+=dt)
    {
      tribol::update(cycle, t, dt);

      op.SetTime(t);
      solver.Step(u, v, t, dt);

      coords += u;
      if (order == 1)
      {
        pmesh->SetVertices(coords);
      }

      ++cycle;
    }

    auto local_max = u.Max();
    max_disp_ = 0.0;
    MPI_Allreduce(&local_max, &max_disp_, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  }
};

TEST_P(MfemCommonPlaneTest, mass_matrix_transfer)
{
  EXPECT_LT(std::abs(max_disp_ - 0.000559868), 1.0e-6);

  MPI_Barrier(MPI_COMM_WORLD);
}

INSTANTIATE_TEST_SUITE_P(tribol, MfemCommonPlaneTest, testing::Values(1));

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

  result = RUN_ALL_TESTS();

  tribol::finalize();
  MPI_Finalize();

  return result;
}