// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Tribol Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (MIT)

// Tribol includes
#include "tribol/types.hpp"
#include "tribol/interface/tribol.hpp"
#include "tribol/utils/TestUtils.hpp"
#include "tribol/utils/Math.hpp"
#include "tribol/common/Parameters.hpp"
#include "tribol/mesh/MethodCouplingData.hpp"
#include "tribol/mesh/CouplingSchemeManager.hpp"
#include "tribol/mesh/CouplingScheme.hpp"
#include "tribol/mesh/InterfacePairs.hpp"
#include "tribol/mesh/MeshData.hpp"
#include "tribol/mesh/MeshManager.hpp"
#include "tribol/geom/GeomUtilities.hpp"

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
 * the computational geometry. This test does not have a specific 
 * check that will make it pass or fail (yet), instead this is 
 * simply used to drive the computational geometry engine 
 * and to interrogate SLIC output printed to screen
 */
class CompGeomTest : public ::testing::Test
{
   
public:

   tribol::TestMesh m_mesh;

protected:

   void SetUp() override
   {
   }

   void TearDown() override
   {
      // call clear() on mesh object to be safe
      this->m_mesh.clear();
   }

protected:

};

TEST_F( CompGeomTest, common_plane_check )
{
   int nMortarElems = 3; 
   int nElemsXM = nMortarElems;
   int nElemsYM = 3;
   int nElemsZM = nMortarElems;

   int nNonmortarElems = 3; 
   int nElemsXS = nNonmortarElems;
   int nElemsYS = 3;
   int nElemsZS = nNonmortarElems;

   int userSpecifiedNumOverlaps = 25;

   // mesh bounding box with 0.1 interpenetration gap
   real x_min1 = 0.;
   real y_min1 = 0.;
   real z_min1 = 0.; 
   real x_max1 = 1.;
   real y_max1 = 1.;
   real z_max1 = 1.05;

   real x_min2 = -0.1;
   real y_min2 = 0.0001;
   real z_min2 = 0.95;
   real x_max2 = 1.1;
   real y_max2 = 0.9999;
   real z_max2 = 2.;

   this->m_mesh.setupContactMeshHex( nElemsXM, nElemsYM, nElemsZM,
                                     x_min1, y_min1, z_min1,
                                     x_max1, y_max1, z_max1,
                                     nElemsXS, nElemsYS, nElemsZS,
                                     x_min2, y_min2, z_min2,
                                     x_max2, y_max2, z_max2,
                                     0., 0. );

   // call tribol setup and update
   tribol::TestControlParameters parameters; // struct does not hold info right now
   parameters.penalty_ratio = false;
   parameters.const_penalty = 1.0;

   int test_mesh_update_err = 
      this->m_mesh.tribolSetupAndUpdate( tribol::COMMON_PLANE, tribol::PENALTY, 
                                         tribol::FRICTIONLESS, true, parameters );
 
   EXPECT_EQ( test_mesh_update_err, 0 );

   tribol::CouplingSchemeManager& couplingSchemeManager = 
         tribol::CouplingSchemeManager::getInstance();
  
   tribol::CouplingScheme* couplingScheme = 
      couplingSchemeManager.getCoupling( 0 );

   EXPECT_EQ( userSpecifiedNumOverlaps, couplingScheme->getNumActivePairs() );
}

TEST_F( CompGeomTest, single_mortar_check )
{
   int nMortarElems = 4; 
   int nElemsXM = nMortarElems;
   int nElemsYM = nMortarElems;
   int nElemsZM = nMortarElems;

   int nNonmortarElems = 5; 
   int nElemsXS = nNonmortarElems;
   int nElemsYS = nNonmortarElems;
   int nElemsZS = nNonmortarElems;

   int userSpecifiedNumOverlaps = 64;

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

   this->m_mesh.setupContactMeshHex( nElemsXM, nElemsYM, nElemsZM,
                                     x_min1, y_min1, z_min1,
                                     x_max1, y_max1, z_max1,
                                     nElemsXS, nElemsYS, nElemsZS,
                                     x_min2, y_min2, z_min2,
                                     x_max2, y_max2, z_max2,
                                     0., 0. );

   // call tribol setup and update
   tribol::TestControlParameters parameters; // struct does not hold info right now

   int test_mesh_update_err = 
      this->m_mesh.tribolSetupAndUpdate( tribol::SINGLE_MORTAR, tribol::LAGRANGE_MULTIPLIER, 
                                         tribol::FRICTIONLESS, false, parameters );

   EXPECT_EQ( test_mesh_update_err, 0 );

   tribol::CouplingSchemeManager& couplingSchemeManager = 
         tribol::CouplingSchemeManager::getInstance();
  
   tribol::CouplingScheme* couplingScheme = 
      couplingSchemeManager.getCoupling( 0 );

   EXPECT_EQ( userSpecifiedNumOverlaps, couplingScheme->getNumActivePairs() );

}

TEST_F( CompGeomTest, poly_area_centroid_1 )
{
   // This test checks the area centroid calculation 
   // vs. the vertex average centroid calculation for a 
   // rectangular quadrilateral. The expectation is that 
   // the results are the same
   int dim = 3;
   int numVerts = 4;
   real x[ dim * numVerts ];

   for (int i=0; i<dim*numVerts; ++i)
   {
      x[i] = 0.;
   }

   // setup some quadrilateral coordinates 
   x[0] = -0.5;
   x[dim*1] =  0.5;
   x[dim*2] =  0.5;
   x[dim*3] = -0.5;

   x[1] = -0.5;
   x[dim*1+1] = -0.5;
   x[dim*2+1] = 0.5;
   x[dim*3+1] = 0.5;

   x[2] = 0.1;
   x[dim*1+2] = 0.1;
   x[dim*2+2] = 0.1; 
   x[dim*3+2] = 0.1;

   real cX_avg, cY_avg, cZ_avg;
   real cX_area, cY_area, cZ_area;

   tribol::VertexAvgCentroid( x, dim, numVerts, cX_avg, cY_avg, cZ_avg );
   tribol::PolyAreaCentroid( x, dim, numVerts, cX_area, cY_area, cZ_area );

   real diff[3] { 0., 0., 0. };

   diff[0] = std::abs(cX_avg - cX_area);
   diff[1] = std::abs(cY_avg - cY_area);
   diff[2] = std::abs(cZ_avg - cZ_area);

   real diff_mag = tribol::magnitude( diff[0], diff[1], diff[2] );

   real tol = 1.e-5;
   EXPECT_LE( diff_mag, tol );
}

TEST_F( CompGeomTest, poly_area_centroid_2 )
{
   // This test checks the area centroid calculation 
   // the centroid calculation for a non-self-intersecting, 
   // closed polygon
   int dim = 3;
   int numVerts = 4;
   real x[ numVerts ];
   real y[ numVerts ];
   real z[ numVerts ];

   for (int i=0; i<numVerts; ++i)
   {
      x[i] = 0.;
      y[i] = 0.;
      z[i] = 0.;
   }

   // setup some quadrilateral coordinates 
   x[0] = -0.515;
   x[1] =  0.54;
   x[2] =  0.65;
   x[3] = -0.524;

   y[0] = -0.5;
   y[1] = -0.5;
   y[2] = 0.5;
   y[3] = 0.5;

   z[0] = 0.1;
   z[1] = 0.1;
   z[2] = 0.1; 
   z[3] = 0.1;

   // create stacked array of coordinates
   real x_bar[ dim*numVerts];
   for (int i=0; i<numVerts; ++i)
   {
      x_bar[dim*i]   = x[i];
      x_bar[dim*i+1] = y[i];
      x_bar[dim*i+2] = z[i];
   }

   real cX_area, cY_area, cZ_area;
   real cX_poly, cY_poly, cZ_poly;

   tribol::PolyAreaCentroid( x_bar, dim, numVerts, cX_area, cY_area, cZ_area );
   tribol::PolyCentroid( x, y, numVerts, cX_poly, cY_poly );

   cZ_poly = z[0];

   real diff[3] { 0., 0., 0. };

   diff[0] = std::abs(cX_poly - cX_area);
   diff[1] = std::abs(cY_poly - cY_area);
   diff[2] = std::abs(cZ_poly - cZ_area);

   real diff_mag = tribol::magnitude( diff[0], diff[1], diff[2] );

   real tol = 1.e-5;
   EXPECT_LE( diff_mag, tol );
}

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger;                // create & initialize logger,
  tribol::SimpleMPIWrapper wrapper(argc, argv);   // initialize and finalize MPI, when applicable

  result = RUN_ALL_TESTS();

  return result;
}
