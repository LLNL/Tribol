// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Tribol Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (MIT)

// Tribol includes
#include "tribol/types.hpp"
#include "tribol/interface/tribol.hpp"
#include "tribol/interface/simple_tribol.hpp"

#include "tribol/utils/TestUtils.hpp"
#include "tribol/utils/Math.hpp"
#include "tribol/common/Parameters.hpp"
#include "tribol/mesh/MethodCouplingData.hpp"
#include "tribol/mesh/CouplingSchemeManager.hpp"
#include "tribol/mesh/CouplingScheme.hpp"
#include "tribol/mesh/MeshData.hpp"
#include "tribol/mesh/MeshManager.hpp"
#include "tribol/physics/Mortar.hpp"
#include "tribol/physics/AlignedMortar.hpp"
#include "tribol/geom/GeomUtilities.hpp"

// Axom includes
#include "axom/slic.hpp"

// MFEM includes
#include "mfem.hpp"

// gtest includes
#include "gtest/gtest.h"

// c++ includes
#include <cmath> // std::abs, std::cos, std::sin
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>

using real = tribol::real;

void computeGapsFromSparseWts( tribol::CouplingScheme const * cs, real * gaps )
{
   tribol::MeshManager& meshManager = tribol::MeshManager::getInstance();
   tribol::parameters_t& parameters = tribol::parameters_t::getInstance();
   tribol::integer const dim = parameters.dimension;

   ////////////////////////////////////////////////////////////////////////
   //
   // Grab pointers to mesh data
   //
   ////////////////////////////////////////////////////////////////////////
   tribol::IndexType const masterId = cs->getMeshId1();
   tribol::IndexType const slaveId = cs->getMeshId2();

   tribol::MeshData& masterMesh = meshManager.GetMeshInstance( masterId );
   tribol::MeshData& slaveMesh = meshManager.GetMeshInstance( slaveId );

   // get mortar weights in CSR format. Note this simple API function 
   // calls tribol::getCSRMatrix() so this API function in the Tribol 
   // namespace is getting tested.
   int *I = nullptr;
   int *J = nullptr;
   real *wts = nullptr;
   int nOffsets = 0;
   int nNonZeros = 0;
   int csr_err = GetSimpleCouplingCSR( &I, &J, &wts, &nOffsets, &nNonZeros );

   EXPECT_EQ( csr_err, 0 );

   if (I == nullptr)
   {
      SLIC_ERROR("Mortar wts test, I is null.");
   }

   // get master node id offset to distinguish master from slave column contributions
   if (masterMesh.m_sortedSurfaceNodeIds == nullptr)
   {
      SLIC_INFO("computeGapsFromSparseWts(): sorting unique master surface node ids.");
      masterMesh.sortSurfaceNodeIds();
   }
   int nodeOffset = masterMesh.m_sortedSurfaceNodeIds[ masterMesh.m_numSurfaceNodes-1 ] + 1;

   ////////////////////////////////////////////////////////////////
   // compute slave gaps to determine active set of contact dofs //
   ////////////////////////////////////////////////////////////////
   // loop over total number of nodes. Note, only the active slave "rows" should 
   // have nonzero contributions
   int numTotalNodes = static_cast<tribol::MortarData*>( cs->getMethodData() )->m_numTotalNodes;
   for (int a=0; a<numTotalNodes; ++a)
   {
      // get slave nodal normal
      real nrml_a[dim];
      nrml_a[0] = slaveMesh.m_node_nX[ a ]; // array is global length; no index out (?)
      nrml_a[1] = slaveMesh.m_node_nY[ a ];
      if (dim == 3 )
      {
         nrml_a[2] = slaveMesh.m_node_nZ[ a ];
      }

      // loop over range of nonzero column entries
      for (int b=I[a]; b<I[a+1]; ++b)
      {
         // get face coordinates for node J[b], i.e. column node id
         real master_xyz[ dim ]; 
         real slave_xyz[ dim ]; 

         for (int i=0; i<dim; ++i)
         {
            master_xyz[i] = 0.;
            slave_xyz[i]  = 0.;
         }

         real n_ab = wts[b];

         if ( J[b] < nodeOffset ) // slave/master  weight
         {
            master_xyz[0] = masterMesh.m_positionX[ J[b] ];
            master_xyz[1] = masterMesh.m_positionY[ J[b] ];
            if ( dim == 3 )
            { 
               master_xyz[2] = masterMesh.m_positionZ[ J[b] ];
            }
 
            gaps[a] += tribol::dotProd( &nrml_a[0], &master_xyz[0], dim ) *
                       n_ab;

         }
         else // slave/slave weight
         {
            slave_xyz[0] = slaveMesh.m_positionX[ J[b] ];
            slave_xyz[1] = slaveMesh.m_positionY[ J[b] ];
            if ( dim == 3 )
            { 
               slave_xyz[2] = slaveMesh.m_positionZ[ J[b] ];
            }
            gaps[a] -= tribol::dotProd( &nrml_a[0], &slave_xyz[0], dim ) * 
                       n_ab;
         } // end if-block

      } // end loop over nonzero columns, I[a]
   } // end loop over matrix rows

} // end ComputeGapsFromSparseWts()

void compareGaps( tribol::CouplingScheme const * cs, real * gaps, const real tol )
{
   tribol::MeshManager& meshManager = tribol::MeshManager::getInstance();
   tribol::IndexType const slaveId = cs->getMeshId2();
   tribol::MeshData& slaveMesh = meshManager.GetMeshInstance( slaveId );

   int numTotalNodes = cs->getNumTotalNodes();

   for (int i=0; i<numTotalNodes; ++i)
   {
      real diff = slaveMesh.m_nodalFields.m_node_gap[i] - gaps[i];
      EXPECT_LE( diff, tol );
   }

}

/*!
 * Test fixture class with some setup necessary to test
 * the MORTAR_WEIGHTS implementation in Tribol
 */
class MortarSparseWtsTest : public ::testing::Test
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

TEST_F( MortarSparseWtsTest, mortar_weights_uniform )
{
   int nMasterElems = 5; 
   int nElemsXM = nMasterElems;
   int nElemsYM = nMasterElems;
   int nElemsZM = nMasterElems;

   int nSlaveElems = 5; 
   int nElemsXS = nSlaveElems;
   int nElemsYS = nSlaveElems;
   int nElemsZS = nSlaveElems;

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
      this->m_mesh.tribolSetupAndUpdate( tribol::MORTAR_WEIGHTS, tribol::NULL_ENFORCEMENT, 
                                         tribol::NULL_MODEL, false, parameters );

   EXPECT_EQ( test_mesh_update_err, 0 );

   tribol::CouplingSchemeManager& couplingSchemeManager = 
         tribol::CouplingSchemeManager::getInstance();
  
   tribol::CouplingScheme* couplingScheme = couplingSchemeManager.getCoupling( 0 );

   // allocate storage for gap computations using sparse mortar weights
   real * gaps = nullptr;
   int size = static_cast<tribol::MortarData*>( couplingScheme->getMethodData() )->m_numTotalNodes;
   gaps = new real[ size ];

   // initialize gap storage
   for (int i=0; i<size; ++i)
   {
      gaps[i] = 0.;
   }

   // use sparse mortar weights to compute gaps
   computeGapsFromSparseWts( couplingScheme, gaps );

   compareGaps( couplingScheme, gaps, 1.E-8 );

   delete [] gaps;
}

TEST_F( MortarSparseWtsTest, simple_api_mortar_weights_uniform )
{
   int nMasterElems = 5; 
   int nElemsXM = nMasterElems;
   int nElemsYM = nMasterElems;
   int nElemsZM = nMasterElems;

   int nSlaveElems = 5; 
   int nElemsXS = nSlaveElems;
   int nElemsYS = nSlaveElems;
   int nElemsZS = nSlaveElems;

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
   int test_mesh_simple_update_err = 
      this->m_mesh.simpleTribolSetupAndUpdate( tribol::MORTAR_WEIGHTS, tribol::NULL_ENFORCEMENT, 
                                               tribol::NULL_MODEL, false, parameters );

   EXPECT_EQ( test_mesh_simple_update_err, 0 );

   tribol::CouplingSchemeManager& couplingSchemeManager = 
         tribol::CouplingSchemeManager::getInstance();
  
   tribol::CouplingScheme* couplingScheme = couplingSchemeManager.getCoupling( 0 );

   // allocate storage for gap computations using sparse mortar weights
   real * gaps = nullptr;
   int size = static_cast<tribol::MortarData*>( couplingScheme->getMethodData() )->m_numTotalNodes;
   gaps = new real[ size ];

   // initialize gap storage
   for (int i=0; i<size; ++i)
   {
      gaps[i] = 0.;
   }

   // use sparse mortar weights to compute gaps
   computeGapsFromSparseWts( couplingScheme, gaps );

   compareGaps( couplingScheme, gaps, 1.E-8 );

   delete [] gaps;
}

TEST_F( MortarSparseWtsTest, mortar_weights_nonuniform_master_fine )
{
   int nMasterElems = 5; 
   int nElemsXM = nMasterElems;
   int nElemsYM = nMasterElems;
   int nElemsZM = nMasterElems;

   int nSlaveElems = 4; 
   int nElemsXS = nSlaveElems;
   int nElemsYS = nSlaveElems;
   int nElemsZS = nSlaveElems;

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
      this->m_mesh.tribolSetupAndUpdate( tribol::MORTAR_WEIGHTS, tribol::NULL_ENFORCEMENT, 
                                         tribol::NULL_MODEL, false, parameters );

   EXPECT_EQ( test_mesh_update_err, 0 );

   tribol::CouplingSchemeManager& couplingSchemeManager = 
         tribol::CouplingSchemeManager::getInstance();
  
   tribol::CouplingScheme* couplingScheme = couplingSchemeManager.getCoupling( 0 );

   // allocate storage for gap computations using sparse mortar weights
   real * gaps = nullptr;
   int size = static_cast<tribol::MortarData*>( couplingScheme->getMethodData() )->m_numTotalNodes;
   gaps = new real[ size ];

   // initialize gap storage
   for (int i=0; i<size; ++i)
   {
      gaps[i] = 0.;
   }

   // use sparse mortar weights to compute gaps
   computeGapsFromSparseWts( couplingScheme, gaps );

   compareGaps( couplingScheme, gaps, 1.E-8 );

   delete [] gaps;
}

TEST_F( MortarSparseWtsTest, mortar_weights_nonuniform_slave_fine )
{
   int nMasterElems = 2; 
   int nElemsXM = nMasterElems;
   int nElemsYM = nMasterElems;
   int nElemsZM = nMasterElems;

   int nSlaveElems = 3; 
   int nElemsXS = nSlaveElems;
   int nElemsYS = nSlaveElems;
   int nElemsZS = nSlaveElems;

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
      this->m_mesh.tribolSetupAndUpdate( tribol::MORTAR_WEIGHTS, tribol::NULL_ENFORCEMENT, 
                                         tribol::NULL_MODEL, false, parameters );

   EXPECT_EQ( test_mesh_update_err, 0 );

   tribol::CouplingSchemeManager& couplingSchemeManager = 
         tribol::CouplingSchemeManager::getInstance();
  
   tribol::CouplingScheme* couplingScheme = couplingSchemeManager.getCoupling( 0 );

   // allocate storage for gap computations using sparse mortar weights
   real * gaps = nullptr;
   int size = static_cast<tribol::MortarData*>( couplingScheme->getMethodData() )->m_numTotalNodes;
   gaps = new real[ size ];

   // initialize gap storage
   for (int i=0; i<size; ++i)
   {
      gaps[i] = 0.;
   }

   // use sparse mortar weights to compute gaps
   computeGapsFromSparseWts( couplingScheme, gaps );

   compareGaps( couplingScheme, gaps, 1.E-8 );

   delete [] gaps;
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
