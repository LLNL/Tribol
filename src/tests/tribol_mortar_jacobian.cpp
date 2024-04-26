// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Tribol Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (MIT)

// Tribol includes
#include "tribol/interface/tribol.hpp"
#include "tribol/common/Parameters.hpp"
#include "tribol/mesh/CouplingScheme.hpp"
#include "tribol/mesh/MethodCouplingData.hpp"
#include "tribol/mesh/MeshData.hpp"
#include "tribol/physics/Mortar.hpp"
#include "tribol/physics/AlignedMortar.hpp"
#include "tribol/geom/GeomUtilities.hpp"
#include "tribol/utils/TestUtils.hpp"

// Axom includes
#include "axom/slic.hpp"

// MFEM includes
#include "mfem.hpp"

// gtest includes
#include "gtest/gtest.h"

// c++ includes
#include <cmath> // std::abs
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>

using RealT = tribol::RealT;

/*!
 * Test fixture class with some setup necessary to test
 * and compute the Jacobian matrix for an implicit 
 * mortar method with Lagrange multipliers
 */
class MortarJacTest : public ::testing::Test
{
   
public:
   int numNodes;
   int numFaces;
   int numNodesPerFace;
   int numOverlapNodes;
   int dim;

   RealT* getXCoords( )
   {
      return x;
   }

   RealT* getYCoords( )
   {
      return y;
   }

   RealT* getZCoords( )
   {
      return z;
   }

   void setupTribol( int * conn1,
                     int * conn2,
                     tribol::ContactMethod method )
   {
      // Note, this assumes that numNodes is the total number of 
      // nodes encompassing the two meshes that will be registered 
      // with tribol, and that the conn1 and conn2 connectivity arrays 
      // reflect a global, contiguous index space

      // grab coordinate data
      RealT * x = this->x;
      RealT * y = this->y;
      RealT * z = this->z;

      // register the mesh with tribol
      int cellType = static_cast<int>(tribol::UNDEFINED_ELEMENT);
      switch (this->numNodesPerFace)
      {
         case 4:
         {
            cellType = (int)(tribol::LINEAR_QUAD);
            break;
         }
         default:
         {
            SLIC_ERROR("setupTribol: number of nodes per face not equal to 4.");
         }
      }

      const int mortarMeshId = 0;
      const int nonmortarMeshId = 1;

      // initialize tribol
      tribol::CommT problem_comm = TRIBOL_COMM_WORLD;
      tribol::initialize( dim, problem_comm );

      // register mesh
      tribol::registerMesh( mortarMeshId, 1, 
                            this->numNodes,
                            conn1, cellType, 
                            x, y, z );
      tribol::registerMesh( nonmortarMeshId, 1, 
                            this->numNodes,
                            conn2, cellType, 
                            x, y, z );

      // register nodal forces. Note, I was getting a seg fault when 
      // registering the same pointer to a single set of force arrays 
      // for both calls to tribol::registerNodalResponse(). As a result,
      // I created two sets of nodal force arrays with their own pointers 
      // to the data that are registered with tribol and there is no longer
      // a seg fault.
      RealT *fx1, *fy1, *fz1;
      RealT *fx2, *fy2, *fz2;

      RealT forceX1[ this->numNodes ];
      RealT forceY1[ this->numNodes ];
      RealT forceZ1[ this->numNodes ];     

      RealT forceX2[ this->numNodes ];
      RealT forceY2[ this->numNodes ];
      RealT forceZ2[ this->numNodes ];     

      fx1 = forceX1; 
      fy1 = forceY1;
      fz1 = forceZ1; 

      fx2 = forceX2;
      fy2 = forceY2;
      fz2 = forceZ2;

      // initialize force arrays
      for (int i=0; i<this->numNodes; ++i)
      {
         fx1[i] = 0.;
         fy1[i] = 0.;
         fz1[i] = 0.;

         fx2[i] = 0.;
         fy2[i] = 0.;
         fz2[i] = 0.;
      }

      tribol::registerNodalResponse( mortarMeshId, fx1, fy1, fz1 );
      tribol::registerNodalResponse( nonmortarMeshId, fx2, fy2, fz2 );

      gaps = tribol::ArrayT<RealT>(this->numNodes, this->numNodes); // length of total mesh to use global connectivity to index
      pressures = tribol::ArrayT<RealT>(this->numNodes, this->numNodes); // length of total mesh to use global connectivity to index

      // initialize gaps and pressures. Initialize all 
      // nonmortar pressures to 1.0
      for (int i=0; i<this->numNodes; ++i)
      {
         gaps[i] = 0.;
         pressures[i] = 1.;
      }

      // register nodal gaps and pressures
      tribol::registerMortarGaps( nonmortarMeshId, gaps.data() );
      tribol::registerMortarPressures( nonmortarMeshId, pressures.data() );

      // register coupling scheme
      const int csIndex = 0;
      tribol::registerCouplingScheme( csIndex,
                                      mortarMeshId,
                                      nonmortarMeshId,
                                      tribol::SURFACE_TO_SURFACE,
                                      tribol::NO_CASE,
                                      method,
                                      tribol::FRICTIONLESS,
                                      tribol::LAGRANGE_MULTIPLIER );

      tribol::setLagrangeMultiplierOptions( csIndex, tribol::ImplicitEvalMode::MORTAR_RESIDUAL_JACOBIAN, 
                                            tribol::SparseMode::MFEM_LINKED_LIST );

      RealT dt = 1.0;
      int tribol_update_err = tribol::update( 1, 1., dt );

      EXPECT_EQ( tribol_update_err, 0 );
 
   }

protected:

   void SetUp() override
   {
      this->numNodes = 8;
      this->numFaces = 2;
      this->numNodesPerFace = 4;
      this->numOverlapNodes = 4;
      this->dim = 3;

      if (this->x == nullptr)
      {
         this->x = new RealT [this->numNodes];
      }
      else
      {
         delete [] this->x;
         this->x = new RealT [this->numNodes];
      }

      if (this->y == nullptr)
      {
         this->y = new RealT [this->numNodes];
      }
      else
      {
         delete [] this->y;
         this->y = new RealT [this->numNodes];
      }

      if (this->z == nullptr)
      {
         this->z = new RealT [this->numNodes];
      }
      else
      {
         delete [] this->z;
         this->z = new RealT [this->numNodes];
      }
  
   }

   void TearDown() override
   {
      if (this->x != nullptr)
      {
         delete [] this->x;
         this->x = nullptr;
      }
      if (this->y != nullptr)
      {
         delete [] this->y;
         this->y = nullptr;
      }
      if (this->z != nullptr)
      {
         delete [] this->z;
         this->z = nullptr;
      }
   }

protected:

   RealT* x {nullptr};
   RealT* y {nullptr};
   RealT* z {nullptr};
   tribol::ArrayT<RealT> gaps;
   tribol::ArrayT<RealT> pressures;

};

TEST_F( MortarJacTest, jac_input_test )
{

   RealT* x = this->getXCoords();
   RealT* y = this->getYCoords();
   RealT* z = this->getZCoords();

   x[0] = -1.;
   x[1] = -1.;
   x[2] =  1.;
   x[3] =  1.;

   y[0] =  1.;
   y[1] = -1.;
   y[2] = -1.;
   y[3] =  1.;

   z[0] = 0.1;
   z[1] = 0.1;
   z[2] = 0.1; 
   z[3] = 0.1;

   x[4] = -1.;
   x[5] =  1.;
   x[6] =  1.;
   x[7] = -1.;

   y[4] =  1.;
   y[5] =  1.;
   y[6] = -1.;
   y[7] = -1.;

   z[4] = 0.;
   z[5] = 0.;
   z[6] = 0.; 
   z[7] = 0.;

   // register a tribol mesh for computing mortar gaps
   int numNodesPerFace = 4;
   int conn1[ numNodesPerFace ];
   int conn2[ numNodesPerFace ];
  
   for (int i=0; i<numNodesPerFace; ++i)
   {
      conn1[i] = i;
      conn2[i] = numNodesPerFace + i;
   }

   this->setupTribol( &conn1[0], &conn2[0], 
                      tribol::ALIGNED_MORTAR );

   // check Jacobian sparse matrix
   mfem::SparseMatrix * jac { nullptr };
   int sparseMatErr = tribol::getJacobianSparseMatrix( &jac, 0 );

   EXPECT_EQ( sparseMatErr, 0 );

   // add each diagonal entry to 1 (i.e. identity matrix) to test adding routine
   for (int i=0; i<2*numNodesPerFace; ++i)
   {
      jac->Add(i, i, 1.);
   }

   // get the number of rows and compare to expected to make sure 
   // sparse matrix initialization is correct
   int my_num_rows = 3 * 2 * numNodesPerFace + 2 * numNodesPerFace;
   int numRows = jac->NumRows();
   EXPECT_EQ( numRows, my_num_rows );

   // get the number of columns and compare to expected
   int numCols = jac->NumCols();
   EXPECT_EQ( numCols, my_num_rows );

   tribol::finalize();

   // delete the jacobian matrix
   delete jac;


}

TEST_F( MortarJacTest, update_jac_test )
{

   RealT* x = this->getXCoords();
   RealT* y = this->getYCoords();
   RealT* z = this->getZCoords();

   x[0] = -1.;
   x[1] = -1.;
   x[2] =  1.;
   x[3] =  1.;

   y[0] =  1.;
   y[1] = -1.;
   y[2] = -1.;
   y[3] =  1.;

   z[0] = 0.1;
   z[1] = 0.1;
   z[2] = 0.1; 
   z[3] = 0.1;

   x[4] = -1.;
   x[5] =  1.;
   x[6] =  1.;
   x[7] = -1.;

   y[4] =  1.;
   y[5] =  1.;
   y[6] = -1.;
   y[7] = -1.;

   z[4] = 0.;
   z[5] = 0.;
   z[6] = 0.; 
   z[7] = 0.;

   // register a tribol mesh for computing mortar gaps
   int numNodesPerFace = 4;
   int conn1[ numNodesPerFace ];
   int conn2[ numNodesPerFace ];
  
   for (int i=0; i<numNodesPerFace; ++i)
   {
      conn1[i] = i;
      conn2[i] = numNodesPerFace + i;
   }

   this->setupTribol( &conn1[0], &conn2[0],
                      tribol::ALIGNED_MORTAR );

   // check Jacobian sparse matrix
   mfem::SparseMatrix * jac { nullptr };
   int sparseMatErr = tribol::getJacobianSparseMatrix( &jac, 0 );

   EXPECT_EQ( sparseMatErr, 0 );

   // make sure sparse matrix is finalized to convert to CSR format
   jac->Finalize();

   int numRows = jac->NumRows();
   int numCols = numRows;

   // convert sparse matrix to dense matrix for output
   mfem::DenseMatrix dJac;
   jac->ToDenseMatrix(dJac);

   std::ofstream matrix;
   matrix.setf(std::ios::scientific);
   matrix.precision(2);
   std::ostringstream suffix_matrix;
   suffix_matrix << "test2" << ".txt";
   matrix.open("matrix_" + suffix_matrix.str());

   for (int i=0; i<numRows; ++i)
   {
      for (int j=0; j<numCols; ++j)
      {
         RealT val = dJac(i,j);
         matrix << val << "  ";
      }
      matrix << "\n";
   }

   matrix.close();

   tribol::finalize();

   // delete the jacobian matrix
   delete jac;

}

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger;                // create & initialize logger,

  result = RUN_ALL_TESTS();

  return result;
}
