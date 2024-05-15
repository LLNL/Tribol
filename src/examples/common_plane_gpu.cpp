// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Tribol Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (MIT)

#include "examples_common.hpp" // for common functionality used in examples

#include "tribol/common/ExecModel.hpp"
#include "tribol/common/Parameters.hpp"
#include "tribol/config.hpp"
#include "tribol/interface/tribol.hpp"

#ifdef TRIBOL_USE_UMPIRE
// Umpire includes
#include "umpire/ResourceManager.hpp"
#endif

// Example command line arguments for running this example. This test creates two rectangular blocks of dimensions (l x w x h). The dimensions are 
// set such that an initial block-intersection exists, triggering the contact interaction. The blocks are discretized per the "block#_res xx yy zz" input 
// arguments (e.g. block1_res 5 3 4, and block2_res 3 2 1), where "xx", "yy", and "zz" are the number of elements in the x, y and z directions. 
// Varying "xx" and "yy" will vary the number of contact overlaps that exist between the two surfaces, and therefore, the amount of contact work.
//
// srun -n1 ./common_plane_ex --block1_res 100 50 10  --block1_min 0 0 0  --block1_max 10 5 1    --block2_res 150 75 10  --block2_min 0 0 0.95  --block2_max 10 5 1.95
// srun -n1 ./common_plane_ex --block1_res 100 50 4   --block1_min 0 0 0. --block1_max 1 1 1.05  --block2_res 150 75 4   --block2_min 0 0 0.95  --block2_max 1 1 2
// srun -n1 ./common_plane_ex --block1_res 4 4 4      --block1_min 0 0 0. --block1_max 1 1 1.05  --block2_res 4 4 4      --block2_min 0 0 0.95  --block2_max 1 1 2

template <tribol::MemorySpace MSPACE, tribol::ExecutionMode EXEC>
int runExample();

/*!
 * \brief Program main.
 *
 * This example runs the common plane + penalty algorithm for two opposing blocks
 *
 * \param [in] argc argument counter
 * \param [in] argv argument vector
 *
 * \return rc return code, a non-zero return code indicates an error.
 */
int main( int argc, char** argv )
{
  ////////////////////////////////
  //                            //
  // SETUP THE EXAMPLE AND MESH //
  //                            //
  ////////////////////////////////

  // initialize
#ifdef TRIBOL_USE_MPI
  MPI_Init( &argc, &argv );
#endif
  tribol::CommT problem_comm = TRIBOL_COMM_WORLD;
  initialize_logger( problem_comm );

#ifdef TRIBOL_USE_UMPIRE
  umpire::ResourceManager::getInstance();  // initialize umpire's ResouceManager
#endif

  auto mem_space = tribol::MemorySpace::Host;

  int err = 0;
  if (mem_space == tribol::MemorySpace::Device)
  {
    err = runExample<tribol::MemorySpace::Device, tribol::ExecutionMode::Cuda>();
  }
  else
  {
    err = runExample<tribol::MemorySpace::Host, tribol::ExecutionMode::Sequential>();
  }

  axom::slic::flushStreams();
  finalize_logger();

#ifdef TRIBOL_USE_MPI
  MPI_Finalize();
#endif

  return err;
}

template <tribol::MemorySpace MSPACE, tribol::ExecutionMode EXEC>
int runExample()
{

  mfem::Device device;
  if (MSPACE == tribol::MemorySpace::Device)
  {
    device.Configure("cuda");
  }
  else
  {
    device.Configure("cpu");
  }
  device.Print();

  std::string mesh_file = TRIBOL_REPO_DIR "/data/two_hex_apart.mesh";
  mfem::Mesh mesh(mesh_file);
  mfem::H1_FECollection fe_coll(1, mesh.SpaceDimension());
  mfem::FiniteElementSpace fe_space(&mesh, &fe_coll, mesh.SpaceDimension());
  mfem::GridFunction coords(&fe_space);
  mesh.SetNodalGridFunction(&coords, false);

  auto coords_ptr = coords.Read();
  auto x_coords_ptr = &coords_ptr[fe_space.DofToVDof(0, 0)];
  auto y_coords_ptr = &coords_ptr[fe_space.DofToVDof(0, 1)];
  auto z_coords_ptr = &coords_ptr[fe_space.DofToVDof(0, 2)];

  // mesh 1 connectivity (build on cpu)
  auto mesh1_bdry_attrib = 4;
  tribol::ArrayT<tribol::IndexT, 2, tribol::MemorySpace::Host> host_conn_1(1, 4);
  for (int be{0}; be < mesh.GetNBE(); ++be)
  {
    if (mesh.GetBdrAttribute(be) == mesh1_bdry_attrib)
    {
      mfem::Array<int> be_dofs(4);//, mfem::MemoryType::Host_UMPIRE);
      fe_space.GetBdrElementDofs(be, be_dofs);
      for (int i{0}; i < 4; ++i)
      {
        host_conn_1(0, i) = be_dofs[i];
      }
    }
  }
  // move to gpu
  tribol::ArrayT<tribol::IndexT, 2, MSPACE> conn_1(host_conn_1);

  // mesh 2 connectivity (build on cpu)
  auto mesh2_bdry_attrib = 5;
  tribol::ArrayT<tribol::IndexT, 2, tribol::MemorySpace::Host> host_conn_2(1, 4);
  for (int be{0}; be < mesh.GetNBE(); ++be)
  {
    if (mesh.GetBdrAttribute(be) == mesh2_bdry_attrib)
    {
      mfem::Array<int> be_dofs(4);//, mfem::MemoryType::Host_UMPIRE);
      fe_space.GetBdrElementDofs(be, be_dofs);
      for (int i{0}; i < 4; ++i)
      {
        host_conn_2(0, i) = be_dofs[i];
      }
    }
  }
  // move to gpu
  tribol::ArrayT<tribol::IndexT, 2, MSPACE> conn_2(host_conn_2);

  constexpr tribol::IndexT mesh1_id = 0;
  tribol::registerMesh(mesh1_id, 1, fe_space.GetNDofs(), conn_1.data(), tribol::LINEAR_QUAD,
    x_coords_ptr, y_coords_ptr, z_coords_ptr, MSPACE);
  constexpr tribol::IndexT mesh2_id = 1;
  tribol::registerMesh(mesh2_id, 1, fe_space.GetNDofs(), conn_2.data(), tribol::LINEAR_QUAD,
    x_coords_ptr, y_coords_ptr, z_coords_ptr, MSPACE);

  constexpr tribol::RealT penalty = 1000.0;
  tribol::setKinematicConstantPenalty(mesh1_id, penalty);
  tribol::setKinematicConstantPenalty(mesh2_id, penalty);

  mfem::GridFunction velocity(&fe_space);
  velocity = 0.0;
  auto velocity_ptr = velocity.Read();
  auto x_velocity_ptr = &velocity_ptr[fe_space.DofToVDof(0, 0)];
  auto y_velocity_ptr = &velocity_ptr[fe_space.DofToVDof(0, 1)];
  auto z_velocity_ptr = &velocity_ptr[fe_space.DofToVDof(0, 2)];
  tribol::registerNodalVelocities(mesh1_id, x_velocity_ptr, y_velocity_ptr, z_velocity_ptr);
  tribol::registerNodalVelocities(mesh2_id, x_velocity_ptr, y_velocity_ptr, z_velocity_ptr);

  mfem::Vector force(fe_space.GetVSize());
  force.UseDevice(true);
  force = 0.0;
  auto force_ptr = force.ReadWrite();
  auto x_force_ptr = &force_ptr[fe_space.DofToVDof(0, 0)];
  auto y_force_ptr = &force_ptr[fe_space.DofToVDof(0, 1)];
  auto z_force_ptr = &force_ptr[fe_space.DofToVDof(0, 2)];
  tribol::registerNodalResponse(mesh1_id, x_force_ptr, y_force_ptr, z_force_ptr);
  tribol::registerNodalResponse(mesh2_id, x_force_ptr, y_force_ptr, z_force_ptr);

  constexpr tribol::IndexT cs_id = 0;
  tribol::registerCouplingScheme(cs_id, mesh1_id, mesh2_id,
    tribol::SURFACE_TO_SURFACE, 
    tribol::NO_CASE, 
    tribol::COMMON_PLANE, 
    tribol::FRICTIONLESS,
    tribol::PENALTY,
    tribol::BINNING_BVH,
    EXEC);

  tribol::setPenaltyOptions(cs_id, tribol::KINEMATIC, tribol::KINEMATIC_CONSTANT);

  int cycle = 1;
  RealT time = 1.0;
  RealT dt = 1.0;
  tribol::update(cycle, time, dt);

  axom::slic::flushStreams();
  finalize_logger();

#ifdef TRIBOL_USE_MPI
  MPI_Finalize();
#endif

  return 0;
}

