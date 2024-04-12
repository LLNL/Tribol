// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Tribol Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (MIT)

// Tribol includes
#include "tribol/types.hpp"
#include "tribol/mesh/MeshData.hpp"
#include "tribol/mesh/MethodCouplingData.hpp"
#include "tribol/integ/Integration.hpp"
#include "tribol/geom/GeomUtilities.hpp"
#include "tribol/integ/FE.hpp"

// Axom includes
#include "axom/slic.hpp"

// gtest includes
#include "gtest/gtest.h"

// c++ includes
#include <cmath> // std::abs

using RealT = tribol::RealT;

/*!
 * Test fixture class with some setup necessary to use the
 * triangular decomposition of a quadrilateral with integration 
 * points specified on each triangle's parent space, forward mapped 
 * to the physical triangle, and then mapped using the inverse 
 * isoparametric mapping in order to obtain (xi,eta) coordinates 
 * on the parent four node quad. These tests compute the area 
 * as calculated by summing integrals of shape functions defined 
 * on the four node quad.
 */
class IsoIntegTest : public ::testing::Test
{
   
public:
   int numNodes;
   int dim;

   RealT* getXCoords() 
   {
      return x;
   }

   RealT* getYCoords() 
   {
      return y;
   }

   RealT* getZCoords() 
   {
      return z;
   }

   bool integrate ( RealT const tol )
   {
      RealT xyz[ this->dim * this->numNodes ];
      RealT* xy = xyz;

      RealT * x = this->x;
      RealT * y = this->y;
      RealT * z = this->z;

      // generate stacked coordinate array
      for (int j=0; j<this->numNodes; ++j)
      {
         for (int k=0; k<this->dim; ++k)
         {
            switch (k)
            {
               case 0:
                  xy[ this->dim * j + k ] = x[ j ];
                  break;
               case 1:
                  xy[ this->dim * j + k ] = y[ j ];
                  break;
               case 2:
                  xy[ this->dim * j + k ] = z[ j ];
                  break;
            } // end switch
         } // end loop over dimension
      } // end loop over nodes

      // instantiate SurfaceContactElem struct. Note, this object is instantiated 
      // using face 1 as face 2, but these faces are not used in this test so this 
      // is ok.
      tribol::SurfaceContactElem elem ( this->dim, xy, xy, xy, 
                                        this->numNodes, this->numNodes, 
                                        0, 1, 0, 0);

      // instantiate integration object
      tribol::IntegPts integ;

      // generate all current configuration integration point coordinates and weights
      tribol::GaussPolyIntTri( elem, integ, 2 );

      // evaluate sum_a (integral_face (phi_a) da) with outer loop over nodes, a, and 
      // inner loop over number of integration points
      RealT areaTest = 0.;
      RealT phi = 0.;

      for (int a=0; a<this->numNodes; ++a)
      {
         for (int ip=0; ip<integ.numIPs; ++ip)
         {
            // perform inverse isoparametric mapping of current configuration 
            // integration point to four node quad parent space
            RealT xp[3] = { integ.xy[dim*ip], integ.xy[dim*ip+1], integ.xy[dim*ip+2] };
            RealT xi[2] = { 0., 0. };
            tribol::InvIso( xp, x, y, z, this->numNodes, xi );
            tribol::LinIsoQuadShapeFunc( xi[0], xi[1], a, phi );
       
            areaTest += integ.wts[ip]*phi;

         }
      }

      RealT area = tribol::Area2DPolygon( x, y, this->numNodes );

      bool convrg = (std::abs(areaTest - area) <= tol) ? true : false;
  
      return convrg;

   }

protected:

   void SetUp() override
   {
      this->numNodes = 4;
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

};


TEST_F( IsoIntegTest, square )
{

   RealT* x = this->getXCoords();
   RealT* y = this->getYCoords();
   RealT* z = this->getZCoords();

   x[0] = -0.5;
   x[1] =  0.5;
   x[2] =  0.5;
   x[3] = -0.5;

   y[0] = -0.5;
   y[1] = -0.5;
   y[2] = 0.5;
   y[3] = 0.5;

   z[0] = 0.1;
   z[1] = 0.1;
   z[2] = 0.1; 
   z[3] = 0.1;

   bool convrg = this->integrate( 1.e-8 );

   EXPECT_EQ( convrg, true );
}

TEST_F( IsoIntegTest, rect )
{

   RealT* x = this->getXCoords();
   RealT* y = this->getYCoords();
   RealT* z = this->getZCoords();

   x[0] = -0.5;
   x[1] =  0.5;
   x[2] =  0.5;
   x[3] = -0.5;

   y[0] = -0.25;
   y[1] = -0.25;
   y[2] = 0.25;
   y[3] = 0.25;

   z[0] = 0.1;
   z[1] = 0.1;
   z[2] = 0.1; 
   z[3] = 0.1;

   bool convrg = this->integrate( 1.e-8 );

   EXPECT_EQ( convrg, true );
}

TEST_F( IsoIntegTest, affine )
{
   RealT* x = this->getXCoords();
   RealT* y = this->getYCoords();
   RealT* z = this->getZCoords();

   x[0] = -0.5;
   x[1] =  0.5;
   x[2] =  0.8;
   x[3] = -0.2;

   y[0] = -0.415;
   y[1] = -0.415;
   y[2] = 0.5;
   y[3] = 0.5;

   z[0] = 0.1;
   z[1] = 0.1;
   z[2] = 0.1; 
   z[3] = 0.1;

   bool convrg = integrate( 1.e-8 );

   EXPECT_EQ( convrg, true );
}

TEST_F( IsoIntegTest, nonaffine )
{
   RealT* x = this->getXCoords();
   RealT* y = this->getYCoords();
   RealT* z = this->getZCoords();

   x[0] = -0.5;
   x[1] =  0.5;
   x[2] =  0.235;
   x[3] = -0.35;

   y[0] = -0.25;
   y[1] = -0.15;
   y[2] = 0.25;
   y[3] = 0.235;

   z[0] = 0.1;
   z[1] = 0.1;
   z[2] = 0.1; 
   z[3] = 0.1;

   // note slightly lower convergence tol for nonaffinely 
   // mapped quad
   bool convrg = integrate( 1.e-5 );

   EXPECT_EQ( convrg, true );
}

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}
