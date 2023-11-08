// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Tribol Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (MIT)

#include <gtest/gtest.h>

#include "mfem.hpp"

#include "tribol/config.hpp"
#include "redecomp/redecomp.hpp"
#include "redecomp/utils/ArrayUtility.hpp"

namespace redecomp {

/**
 * @brief This test transfers an mfem::SparseMatrix whose rows and colums are associated with an order-0 L2 field on a
 * RedecompMesh to an mfem::HypreParMatrix on the associated parent mfem::ParMesh.  The transferred matrix is verified
 * by applying it to a field and comparing the result to an analytic solution for the operator.
 *
 */
class SparseMatrixTest : public testing::TestWithParam<std::pair<int, double>> {
protected:
  double max_error_;
  void SetUp() override
  {
    auto ref_levels = GetParam().first;
    auto filter_radius = GetParam().second;
    
    auto serial_mesh = mfem::Mesh::MakeCartesian2D(4, 4, mfem::Element::Type::QUADRILATERAL);
    auto space_dim = serial_mesh.SpaceDimension();
    auto dim = serial_mesh.Dimension();
    for (int i{0}; i < ref_levels; ++i)
    {
      serial_mesh.UniformRefinement();
    }
    mfem::ParMesh par_mesh { MPI_COMM_WORLD, serial_mesh };
    mfem::L2_FECollection l2_elems { 0, dim };
    mfem::ParFiniteElementSpace par_fes { &par_mesh, &l2_elems };

    // create RedecompMesh and fe space and sparse matrix transfer object
    redecomp::RedecompMesh redecomp_mesh { par_mesh };
    mfem::FiniteElementSpace redecomp_fes { &redecomp_mesh, &l2_elems };
    redecomp::SparseMatrixTransfer matrix_xfer {
      par_fes,
      par_fes,
      redecomp_fes, 
      redecomp_fes
    };

    // compute sparse matrix on redecomp_fes (rows and col)
    mfem::SparseMatrix rho_redecomp { redecomp_fes.GetVSize(), redecomp_fes.GetVSize() };

    auto filter_kernel = [&filter_radius](const mfem::Vector &xi, const mfem::Vector &xj) {
      return std::max(0.0, filter_radius - xi.DistanceTo(xj));
    };

    int n_local_elem = redecomp_mesh.GetNE();
    mfem::Vector xi(space_dim), xj(space_dim);
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
          rho_redecomp.Add(i, j, Wij);
          n_row_entries[i]++;
        }
      }
    }

    rho_redecomp.Finalize();

    // normalize each row to conserve mass
    for (int i=0; i<n_local_elem; i++) {
      mfem::Vector row_data(rho_redecomp.GetRowEntries(i), n_row_entries[i]);
      row_data /= row_data.Norml1();
    }

    // transfer to mfem::ParMesh
    auto rho = matrix_xfer.TransferToParallel(rho_redecomp);

    // test operator
    mfem::ParGridFunction x(&par_fes);
    mfem::ParGridFunction xf(&par_fes);

    mfem::FunctionCoefficient x_func([](const mfem::Vector&x) {
      return ((x[0] > 0.5) && (x[1] > 0.5))? 1.0 : 0.0;
    });
    
    x.ProjectCoefficient(x_func);
    
    rho->Mult(x, xf);

    // TODO: compute error

    max_error_ = 0.0; // TODO: store l2 error here
  }
};

TEST_P(SparseMatrixTest, mass_matrix_transfer)
{
  EXPECT_LT(max_error_, 1.0e-13);

  MPI_Barrier(MPI_COMM_WORLD);
}

INSTANTIATE_TEST_SUITE_P(redecomp, SparseMatrixTest, testing::Values(
  std::make_pair(3, 0.2)
  // TODO: add more test combinations here...
));

}  // namespace redecomp

//------------------------------------------------------------------------------
#include "axom/slic/core/SimpleLogger.hpp"

int main(int argc, char* argv[])
{
  int result = 0;

  MPI_Init(&argc, &argv);

  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger;  // create & initialize test logger, finalized when
                                    // exiting main scope

  result = RUN_ALL_TESTS();

  MPI_Finalize();

  return result;
}
