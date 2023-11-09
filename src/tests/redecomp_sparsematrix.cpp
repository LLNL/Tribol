// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Tribol Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (MIT)

#include <gtest/gtest.h>

#include "mfem.hpp"

#include "redecomp/RedecompMesh.hpp"
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
class SparseMatrixTest : public testing::TestWithParam<std::tuple<int, int, double>> {
protected:
  double max_error_;
  void SetUp() override
  {
    auto dim = std::get<0>(GetParam());
    auto ref_levels = std::get<1>(GetParam());
    auto filter_radius = std::get<2>(GetParam());
    
    double side_length = 1.0;
    int n_elem_per_dim = 4;
    mfem::Mesh serial_mesh;
    if (dim == 2) {
      serial_mesh = mfem::Mesh::MakeCartesian2D(n_elem_per_dim, n_elem_per_dim, 
        mfem::Element::Type::QUADRILATERAL, false, side_length, side_length);
    }
    else {
      serial_mesh = mfem::Mesh::MakeCartesian3D(n_elem_per_dim, n_elem_per_dim, n_elem_per_dim, 
        mfem::Element::Type::HEXAHEDRON, side_length, side_length, side_length);
    }
    // refine mesh
    for (int i{0}; i < ref_levels; ++i){
      serial_mesh.UniformRefinement();
    }
    n_elem_per_dim *= std::pow(2, ref_levels);
    double dx = side_length / n_elem_per_dim;

    mfem::ParMesh par_mesh { MPI_COMM_WORLD, serial_mesh };
    mfem::L2_FECollection l2_elems { 0, dim };
    mfem::ParFiniteElementSpace par_fes { &par_mesh, &l2_elems };

    // create RedecompMesh and fe space and sparse matrix transfer object
    redecomp::RedecompMesh redecomp_mesh { par_mesh, redecomp::RedecompMesh::RCB, filter_radius };
    mfem::FiniteElementSpace redecomp_fes { &redecomp_mesh, &l2_elems };
    redecomp::SparseMatrixTransfer matrix_xfer {
      par_fes,
      par_fes,
      redecomp_fes, 
      redecomp_fes
    };

    // compute sparse matrix on redecomp_fes (rows and col)
    mfem::SparseMatrix W_redecomp { redecomp_fes.GetVSize(), redecomp_fes.GetVSize() };

    auto filter_kernel = [&filter_radius](const mfem::Vector &xi, const mfem::Vector &xj) {
      return std::max(0.0, filter_radius - xi.DistanceTo(xj));
    };

    int n_local_elem = redecomp_mesh.GetNE();
    mfem::Vector xi(dim), xj(dim);
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
          W_redecomp.Add(i, j, Wij);
          n_row_entries[i]++;
        }
      }
    }

    W_redecomp.Finalize();

    // normalize each row to conserve mass
    for (int i=0; i<n_local_elem; i++) {
      mfem::Vector row_data(W_redecomp.GetRowEntries(i), n_row_entries[i]);
      row_data /= row_data.Norml1();
    }

    // transfer to mfem::ParMesh
    auto W = matrix_xfer.TransferToParallel(W_redecomp);

    // test operator
    auto x_function = [&](const mfem::Vector&x) {
      double pi = 3.1415;
      double f = 0.0;
      for (int i=0; i<dim; i++) {
        f += std::sin(2.*pi*x[i]/side_length);
      }
      return 0.5 + 0.5*f/dim;
    };

    auto xf_function = [&](const mfem::Vector&x) {
      mfem::Vector xj(dim);
      double value = 0.0;
      double denom = 0.0;
      if (dim == 2) {
        for (int i=0; i<n_elem_per_dim; i++) {
          xj(0) = dx/2. + i*dx;
          for (int j=0; j<n_elem_per_dim; j++) {
            xj(1) = dx/2. + j*dx;
            value += filter_kernel(x, xj)*x_function(xj);
            denom += filter_kernel(x, xj);
          }
        }
      }
      else {
        for (int i=0; i<n_elem_per_dim; i++) {
          xj(0) = dx/2. + i*dx;
          for (int j=0; j<n_elem_per_dim; j++) {
            xj(1) = dx/2. + j*dx;
              for (int k=0; k<n_elem_per_dim; k++) {
                xj(2) = dx/2. + k*dx;
                value += filter_kernel(x, xj)*x_function(xj);
                denom += filter_kernel(x, xj);
            }
          }
        }
      }
      return value / denom;
    };

    mfem::FunctionCoefficient xCoef(x_function);
    mfem::FunctionCoefficient xfCoef(xf_function);
    
    // compare filtered filed with analytical solution
    mfem::ParGridFunction x(&par_fes);
    mfem::ParGridFunction xf_filt(&par_fes);
    mfem::ParGridFunction xf_func(&par_fes);

    x.ProjectCoefficient(xCoef);
    W->Mult(x, xf_filt);
    xf_func.ProjectCoefficient(xfCoef);

    // compute error
    mfem::ParGridFunction error  = xf_filt;
                          error -= xf_func;
    max_error_ = mfem::InnerProduct(par_mesh.GetComm(), error, error);

    // debug output
    // mfem::VisItDataCollection redecomp_dc { "redecomp_rank" + std::to_string(par_mesh.GetMyRank()), &redecomp_mesh };
    // redecomp_dc.Save();

    // mfem::VisItDataCollection visit_dc("test_look", &par_mesh);
    // visit_dc.RegisterField("x", &x);
    // visit_dc.RegisterField("xf_filt", &xf_filt);
    // visit_dc.RegisterField("xf_func", &xf_func);
    // visit_dc.RegisterField("error", &error);
    // visit_dc.Save();

    // std::cout << "filt" << std::endl;
    // xf_filt.Print();
    // std::cout << "func" << std::endl;
    // xf_func.Print();
    // std::cout << "error" << std::endl;
    // error.Print();
  }
};

TEST_P(SparseMatrixTest, mass_matrix_transfer)
{
  EXPECT_LT(max_error_, 1.0e-13);

  MPI_Barrier(MPI_COMM_WORLD);
}

INSTANTIATE_TEST_SUITE_P(redecomp, SparseMatrixTest, testing::Values(
  std::make_tuple(2, 0, 0.1),
  std::make_tuple(2, 0, 0.4),

  std::make_tuple(2, 1, 0.1),
  std::make_tuple(2, 1, 0.4),

  std::make_tuple(2, 2, 0.1),
  std::make_tuple(2, 2, 0.4),
  std::make_tuple(2, 3, 0.1),
  std::make_tuple(2, 3, 0.4),

  std::make_tuple(3, 0, 0.1),
  std::make_tuple(3, 0, 0.4),
  std::make_tuple(3, 1, 0.1),
  std::make_tuple(3, 1, 0.4),
  std::make_tuple(3, 2, 0.1),
  std::make_tuple(3, 2, 0.4)
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
