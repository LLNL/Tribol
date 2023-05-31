// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Tribol Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (MIT)

// Tribol includes
#include "tribol/types.hpp"
#include "tribol/utils/TestUtils.hpp"
#include "tribol/utils/Math.hpp"
#include "tribol/common/Parameters.hpp"

// Axom includes
#include "axom/slic.hpp"

// gtest includes
#include "gtest/gtest.h"

// c++ includes
#include <cmath> // std::abs, std::cos, std::sin
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>

using real = tribol::real;

/*!
 * Test fixture class with some setup necessary to test 
 * the TestMesh tet mesh feature 
 */
class TetMeshTest : public ::testing::Test
{
   
public:

   tribol::TestMesh m_mesh;

protected:

   void SetUp() override
   {
   }

   void TearDown() override
   {
      this->m_mesh.clear();
   }

protected:

};

TEST_F( TetMeshTest, build_and_check_tet_mesh )
{
   this->m_mesh.mortarMeshId = 0;
   this->m_mesh.nonmortarMeshId = 1;

   int nMortarElems = 4;
   int nElemsXM = nMortarElems;
   int nElemsYM = nMortarElems;
   int nElemsZM = nMortarElems;

   int nNonmortarElems = 3; 
   int nElemsXS = nNonmortarElems;
   int nElemsYS = nNonmortarElems;
   int nElemsZS = nNonmortarElems;

   // mesh bounding box with 0.1 interpenetration gap
   real x_min1 = 0.;
   real y_min1 = 0.;
   real z_min1 = 0.; 
   real x_max1 = 1.;
   real y_max1 = 1.;
   real z_max1 = 1.05;

   real x_min2 = 0.;
   real y_min2 = 0.;
   real z_min2 = 0.95;
   real x_max2 = 1.;
   real y_max2 = 1.;
   real z_max2 = 2.;

   this->m_mesh.setupContactMeshTet( nElemsXM, nElemsYM, nElemsZM,
                                     x_min1, y_min1, z_min1,
                                     x_max1, y_max1, z_max1,
                                     nElemsXS, nElemsYS, nElemsZS,
                                     x_min2, y_min2, z_min2,
                                     x_max2, y_max2, z_max2,
                                     5., -5. );
    
   // sanity checks for this specific mesh
   EXPECT_EQ( this->m_mesh.mesh_constructed, true );
   EXPECT_EQ( this->m_mesh.cellType, (int)(tribol::LINEAR_TRIANGLE) );
   EXPECT_EQ( this->m_mesh.numNodesPerFace, 3 );
   EXPECT_EQ( this->m_mesh.numNodesPerElement, 4 );
   EXPECT_EQ( this->m_mesh.numMortarElements, 6 * nMortarElems*nMortarElems*nMortarElems );
   EXPECT_EQ( this->m_mesh.numNonmortarElements, 6 * nNonmortarElems*nNonmortarElems*nNonmortarElems );
   EXPECT_EQ( this->m_mesh.numTotalElements, this->m_mesh.numMortarElements + this->m_mesh.numNonmortarElements );
   EXPECT_EQ( this->m_mesh.numMortarFaces, 2 * nMortarElems*nMortarElems );
   EXPECT_EQ( this->m_mesh.numNonmortarFaces, 2 * nNonmortarElems*nNonmortarElems );
   EXPECT_EQ( this->m_mesh.numMortarNodes, (nMortarElems+1)*(nMortarElems+1)*(nMortarElems+1) );
   EXPECT_EQ( this->m_mesh.numNonmortarNodes, (nNonmortarElems+1)*(nNonmortarElems+1)*(nNonmortarElems+1) );

   //this->m_mesh.testMeshToVtk( "", 1, 1 );

}

//TEST_F( TetMeshTest, build_and_check_mfem_tet_mesh )
//{
//   this->m_mesh.mortarMeshId = 0;
//   this->m_mesh.nonmortarMeshId = 1;
//
//   int nMortarElems = 4; 
//   int nElemsXM = nMortarElems;
//   int nElemsYM = nMortarElems;
//   int nElemsZM = nMortarElems;
//
//   int nNonmortarElems = 5; 
//   int nElemsXS = nNonmortarElems;
//   int nElemsYS = nNonmortarElems;
//   int nElemsZS = nNonmortarElems;
//
//   // mesh bounding box with 0.1 interpenetration gap
//   real x_min1 = 0.;
//   real y_min1 = 0.;
//   real z_min1 = 0.; 
//   real x_max1 = 1.;
//   real y_max1 = 1.;
//   real z_max1 = 1.05;
//
//   real x_min2 = 0.;
//   real y_min2 = 0.;
//   real z_min2 = 0.95;
//   real x_max2 = 1.;
//   real y_max2 = 1.;
//   real z_max2 = 2.;
//
//   this->m_mesh.setupContactMeshTet( nElemsXM, nElemsYM, nElemsZM,
//                                     x_min1, y_min1, z_min1,
//                                     x_max1, y_max1, z_max1,
//                                     nElemsXS, nElemsYS, nElemsZS,
//                                     x_min2, y_min2, z_min2,
//                                     x_max2, y_max2, z_max2,
//                                     0., 0. );
//
//   this->m_mesh.setupMfemMesh();
//
//}

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger;
  result = RUN_ALL_TESTS();

  return result;
}
