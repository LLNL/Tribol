// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Tribol Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (MIT)

// Tribol includes
#include "tribol/types.hpp"
#include "tribol/mesh/MeshData.hpp"
#include "tribol/mesh/MethodCouplingData.hpp"
#include "tribol/physics/Mortar.hpp"
#include "tribol/geom/GeomUtilities.hpp"

// Axom includes
#include "axom/slic.hpp"

// gtest includes
#include "gtest/gtest.h"

// c++ includes
#include <cmath> // std::abs

using real = tribol::real;

/*!
 * Test fixture class with some setup necessary to compute 
 * the slave/slave and master/slave mortar projection weights 
 * in order to test a simple projection between two 
 * mis-aligned faces.
 */
class MortarWeightTest : public ::testing::Test
{
   
public:
   int numNodes;
   int numNodesPerFace;
   int numOverlapNodes;
   int dim;

   real* getXCoords( int id ) 
   {
      if (id == 0)
      {
         return x1;
      }
      else 
      {
         return x2;
      }
   }

   real* getYCoords( int id ) 
   {
      if (id == 0)
      {
         return y1;
      }
      else 
      {
         return y2;
      }
   }

   real* getZCoords( int id ) 
   {
      if (id == 0)
      {
         return z1;
      }
      else
      {
         return z2;
      }
   }

   real* getXOverlapCoords()
   {
      return xOverlap;
   }

   real* getYOverlapCoords()
   {
      return yOverlap;
   }

   real* getZOverlapCoords()
   {
      return zOverlap;
   }

   void checkMortarWts ( real pTestSlave[4], real pTestMaster[4])
   {
      if (this->numNodesPerFace != 4)
      {
         SLIC_ERROR("checkMortarWts: number of nodes per face not equal to 4.");
      }

      real xyz1[ this->dim * this->numNodesPerFace ];
      real xyz2[ this->dim * this->numNodesPerFace ];
      real xyzOverlap[ this->dim * this->numOverlapNodes ];
      real* xy1 = xyz1;
      real* xy2 = xyz2;
      real* xyOverlap = xyzOverlap;

      real * x1 = this->x1;
      real * y1 = this->y1;
      real * z1 = this->z1;
      real * x2 = this->x2;
      real * y2 = this->y2;
      real * z2 = this->z2;
      real * xo = this->xOverlap;
      real * yo = this->yOverlap;
      real * zo = this->zOverlap;

      // generate stacked coordinate array
      for (int j=0; j<this->numNodesPerFace; ++j)
      {
         for (int k=0; k<this->dim; ++k)
         {
            switch (k)
            {
               case 0:
                  xy1[ this->dim * j + k ] = x1[ j ];
                  xy2[ this->dim * j + k ] = x2[ j ];
                  xyOverlap[ this->dim * j + k ] = xo[ j ];
                  break;
               case 1:
                  xy1[ this->dim * j + k ] = y1[ j ];
                  xy2[ this->dim * j + k ] = y2[ j ];
                  xyOverlap[ this->dim * j + k ] = yo[ j ];
                  break;
               case 2:
                  xy1[ this->dim * j + k ] = z1[ j ];
                  xy2[ this->dim * j + k ] = z2[ j ];
                  xyOverlap[ this->dim * j + k ] = zo[ j ];
                  break;
            } // end switch
         } // end loop over dimension
      } // end loop over nodes

      // instantiate SurfaceContactElem struct. Note, this object is instantiated 
      // using face 1, face 2, and the set overlap polygon. Note, the mesh ids are set 
      // equal to 0, and the face ids are 0 and 1, respectively.
      tribol::SurfaceContactElem elem ( this->dim, xy1, xy2, xyOverlap, 
                                        this->numNodesPerFace, this->numOverlapNodes, 
                                        0, 0, 0, 1);

      // compute the mortar weights to be stored on 
      // the surface contact element struct.
      tribol::ComputeMortarWeights( elem );

      // test the projection on an arbitrary vector
      real p[4] = {1., 1., 1., 1.};
      for (int i=0; i<this->numNodesPerFace; ++i)
      {
         for (int j=0; j<this->numNodesPerFace; ++j)
         {
            int slaveId = this->numNodesPerFace*i+j;
            int masterId = this->numNodesPerFace * this->numNodesPerFace + 
                           this->numNodesPerFace * i + j;
            pTestSlave[i] += elem.mortarWts[ slaveId ] * p[j];  
            pTestMaster[i] += elem.mortarWts[ masterId ] * p[j];
         }
      }

   }

protected:

   void SetUp() override
   {
      this->numNodes = 8;
      this->numNodesPerFace = 4;
      this->numOverlapNodes = 4;
      this->dim = 3;

      if (this->x1 == nullptr)
      {
         this->x1 = new real [this->numNodes];
      }
      else
      {
         delete [] this->x1;
         this->x1 = new real [this->numNodes];
      }

      if (this->x2 == nullptr)
      {
         this->x2 = new real [this->numNodes];
      }
      else
      {
         delete [] this->x2;
         this->x2 = new real [this->numNodes];
      }

      if (this->y1 == nullptr)
      {
         this->y1 = new real [this->numNodes];
      }
      else
      {
         delete [] this->y1;
         this->y1 = new real [this->numNodes];
      }

      if (this->y2 == nullptr)
      {
         this->y2 = new real [this->numNodes];
      }
      else
      {
         delete [] this->y2;
         this->y2 = new real [this->numNodes];
      }

      if (this->z1 == nullptr)
      {
         this->z1 = new real [this->numNodes];
      }
      else
      {
         delete [] this->z1;
         this->z1 = new real [this->numNodes];
      }

      if (this->z2 == nullptr)
      {
         this->z2 = new real [this->numNodes];
      }
      else
      {
         delete [] this->z2;
         this->z2 = new real [this->numNodes];
      }
  
      if (this->xOverlap == nullptr)
      {
         this->xOverlap = new real [this->numOverlapNodes];
      }
      else
      {
         delete [] this->xOverlap;
         this->xOverlap = new real [this->numOverlapNodes];
      }

      if (this->yOverlap == nullptr)
      {
         this->yOverlap = new real [this->numOverlapNodes];
      }
      else
      {
         delete [] this->yOverlap;
         this->yOverlap = new real [this->numOverlapNodes];
      }
      if (this->zOverlap == nullptr)
      {
         this->zOverlap = new real [this->numOverlapNodes];
      }
      else
      {
         delete [] this->zOverlap;
         this->zOverlap = new real [this->numOverlapNodes];
      }
   }

   void TearDown() override
   {
      if (this->x1 != nullptr)
      {
         delete [] this->x1;
         this->x1 = nullptr;
      }
      if (this->x2 != nullptr)
      {
         delete [] this->x2;
         this->x2 = nullptr;
      }
      if (this->y1 != nullptr)
      {
         delete [] this->y1;
         this->y1 = nullptr;
      }
      if (this->y2 != nullptr)
      {
         delete [] this->y2;
         this->y2 = nullptr;
      }
      if (this->z1 != nullptr)
      {
         delete [] this->z1;
         this->z1 = nullptr;
      }
      if (this->z2 != nullptr)
      {
         delete [] this->z2;
         this->z2 = nullptr;
      }
      if (this->xOverlap != nullptr)
      {
         delete [] this->xOverlap;
         this->xOverlap = nullptr;
      }
      if (this->yOverlap != nullptr)
      {
         delete [] this->yOverlap;
         this->yOverlap = nullptr;
      }
      if (this->zOverlap != nullptr)
      {
         delete [] this->zOverlap;
         this->zOverlap = nullptr;
      }
   }

protected:

   real* x1 {nullptr};
   real* y1 {nullptr};
   real* z1 {nullptr};

   real* x2 {nullptr};
   real* y2 {nullptr};
   real* z2 {nullptr};

   real* xOverlap {nullptr};
   real* yOverlap {nullptr};
   real* zOverlap {nullptr};

};

TEST_F( MortarWeightTest, simple_projection )
{

   real* x1 = this->getXCoords(0);
   real* y1 = this->getYCoords(0);
   real* z1 = this->getZCoords(0);

   real* x2 = this->getXCoords(1);
   real* y2 = this->getYCoords(1);
   real* z2 = this->getZCoords(1);

   real* xOvrlp = this->getXOverlapCoords();
   real* yOvrlp = this->getYOverlapCoords();
   real* zOvrlp = this->getZOverlapCoords();

   x1[0] = -1.;
   x1[1] = -1.;
   x1[2] =  1.;
   x1[3] =  1.;

   y1[0] =  1.;
   y1[1] = -1.;
   y1[2] = -1.;
   y1[3] =  1.;

   z1[0] = 0.1;
   z1[1] = 0.1;
   z1[2] = 0.1; 
   z1[3] = 0.1;

   x2[0] =  0.;
   x2[1] =  2.;
   x2[2] =  2.;
   x2[3] =  0.;

   y2[0] =  0.;
   y2[1] =  0.;
   y2[2] = -2.;
   y2[3] = -2.;

   z2[0] = 0.1;
   z2[1] = 0.1;
   z2[2] = 0.1; 
   z2[3] = 0.1;

   xOvrlp[0] = 0.;
   xOvrlp[1] = 0.;
   xOvrlp[2] = 1.;
   xOvrlp[3] = 1.;

   yOvrlp[0] = 0.;
   yOvrlp[1] = -1.;
   yOvrlp[2] = -1.;
   yOvrlp[3] = 0.;

   zOvrlp[0] = 0.1;
   zOvrlp[1] = 0.1;
   zOvrlp[2] = 0.1;
   zOvrlp[3] = 0.1;

   real pSlave[4] = {0., 0., 0., 0.};
   real pMaster[4] = {0., 0., 0., 0};
   this->checkMortarWts( pSlave, pMaster );

   // hard-code diffs for each element in pSlave and pMaster. 
   // Note, these hard coded values DEPEND on the ordering of 
   // the nodes in this function (above).
   real diffSlave1 = std::abs(0.5625 - pSlave[0]);
   real diffSlave2 = std::abs(0.1875 - pSlave[1]);
   real diffSlave3 = std::abs(0.0625 - pSlave[2]);
   real diffSlave4 = std::abs(0.1875 - pSlave[3]);

   real diffMaster1 = std::abs(0.0625 - pMaster[0]);
   real diffMaster2 = std::abs(0.1875 - pMaster[1]);
   real diffMaster3 = std::abs(0.5625 - pMaster[2]);
   real diffMaster4 = std::abs(0.1875 - pMaster[3]);

   real tol = 1.e-8;
   EXPECT_LE( diffSlave1, tol );
   EXPECT_LE( diffSlave2, tol );
   EXPECT_LE( diffSlave3, tol );
   EXPECT_LE( diffSlave4, tol );

   EXPECT_LE( diffMaster1, tol );
   EXPECT_LE( diffMaster2, tol );
   EXPECT_LE( diffMaster3, tol );
   EXPECT_LE( diffMaster4, tol );

}

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}
