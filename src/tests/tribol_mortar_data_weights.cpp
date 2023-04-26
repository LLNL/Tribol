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
#include "axom/core.hpp"
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
namespace axom_fs = axom::utilities::filesystem;

void TestMortarWeights( tribol::CouplingScheme const * cs, double exact_area, double tol )
{
   ////////////////////////////////////////////////////////////////////////
   //
   // Grab pointers to mesh data
   //
   ////////////////////////////////////////////////////////////////////////
   tribol::MeshManager& meshManager = tribol::MeshManager::getInstance();
   tribol::IndexType const masterId = cs->getMeshId1();
   //tribol::IndexType const slaveId = cs->getMeshId2();

   tribol::MeshData& masterMesh = meshManager.GetMeshInstance( masterId );
   //tribol::MeshData& slaveMesh = meshManager.GetMeshInstance( slaveId );

   // get CSR weights data
   int *I = nullptr;
   int *J = nullptr;
   double *wts = nullptr;
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

   double area = 0.;
   int numTotalNodes = static_cast<tribol::MortarData*>( cs->getMethodData() )->m_numTotalNodes;
   for (int a=0; a<numTotalNodes; ++a)
   {
      // loop over range of nonzero column entries
      for (int b=I[a]; b<I[a+1]; ++b)
      {
         area += wts[b];
          
         if ( J[b] < nodeOffset ) // slave/master  weight
         {
//            SLIC_INFO("slave/master weight for slave node, " << a << " and master node, " << J[b] << ".");
         }
         else // slave/slave weight
         {
//            SLIC_INFO("slave/slave weight for slave node, " << a << " and master node, " << J[b] << ".");
         } // end if-block

      } // end loop over nonzero columns, I[a]
   } // end loop over matrix rows

   area /= 2.;

   SLIC_INFO("area: " << area << ".");

   double diff = std::abs( area - exact_area );
   EXPECT_LE( diff, tol );

} // end TestMortarWeights()
   

/*!
 * Test fixture class with some setup necessary to test
 * the MORTAR_WEIGHTS implementation in Tribol
 */
class MortarSparseWtsTest : public ::testing::Test
{
   
public:
   mfem::Vector v_xm;
   mfem::Vector v_ym; 
   mfem::Vector v_zm; 
   mfem::Vector v_xs; 
   mfem::Vector v_ys; 
   mfem::Vector v_zs;
   mfem::Array<int> v_ixm;
   mfem::Array<int> v_ixs;

   int lengthMasterConn;
   int lengthSlaveConn;
   int lengthMasterNodes;
   int lengthSlaveNodes;
   int numMasterCells;
   int numSlaveCells;

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

TEST_F( MortarSparseWtsTest, mortar_sphere )
{

   // read data sets for mesh
   std::string ixm_file = axom_fs::joinPath(TRIBOL_DATA_DIR,"mortar_sphere/ixm.txt");
   std::string ixs_file = axom_fs::joinPath(TRIBOL_DATA_DIR,"mortar_sphere/ixs.txt");
   std::string xm_file =  axom_fs::joinPath(TRIBOL_DATA_DIR,"mortar_sphere/xm.txt");
   std::string ym_file =  axom_fs::joinPath(TRIBOL_DATA_DIR,"mortar_sphere/ym.txt");
   std::string zm_file =  axom_fs::joinPath(TRIBOL_DATA_DIR,"mortar_sphere/zm.txt");
   std::string xs_file =  axom_fs::joinPath(TRIBOL_DATA_DIR,"mortar_sphere/xs.txt");
   std::string ys_file =  axom_fs::joinPath(TRIBOL_DATA_DIR,"mortar_sphere/ys.txt");
   std::string zs_file =  axom_fs::joinPath(TRIBOL_DATA_DIR,"mortar_sphere/zs.txt");

   this->lengthMasterConn  = 192;
   this->lengthSlaveConn   = 768;
   this->lengthMasterNodes = 834;
   this->lengthSlaveNodes  = 834;
   this->numMasterCells = 48;
   this->numSlaveCells = 192;

   std::ifstream i_ixm( ixm_file ); 
   std::ifstream i_ixs( ixs_file );
   std::ifstream i_xm( xm_file );
   std::ifstream i_ym( ym_file );
   std::ifstream i_zm( zm_file );
   std::ifstream i_xs( xs_file );
   std::ifstream i_ys( ys_file );
   std::ifstream i_zs( zs_file );

   this->v_ixm.SetSize( this->lengthMasterConn );
   this->v_ixs.SetSize( this->lengthSlaveConn );

   this->v_ixm.Load( i_ixm, 1 );
   this->v_ixs.Load( i_ixs, 1 );
   this->v_xm.Load( i_xm, this->lengthMasterNodes );
   this->v_ym.Load( i_ym, this->lengthMasterNodes );
   this->v_zm.Load( i_zm, this->lengthMasterNodes );
   this->v_xs.Load( i_xs, this->lengthSlaveNodes );
   this->v_ys.Load( i_ys, this->lengthSlaveNodes );
   this->v_zs.Load( i_zs, this->lengthSlaveNodes );

   i_ixm.close();
   i_ixs.close();
   i_xm.close();
   i_ym.close();
   i_zm.close();
   i_xs.close();
   i_ys.close();
   i_zs.close();

   SLIC_INFO("After loading mesh data and constructing mfem vectors.");

   // get pointers to mfem vector data
   int* ixm_data   = this->v_ixm.GetData();
   int* ixs_data   = this->v_ixs.GetData();
   double* xm_data = this->v_xm.GetData();
   double* ym_data = this->v_ym.GetData();
   double* zm_data = this->v_zm.GetData();
   double* xs_data = this->v_xs.GetData();
   double* ys_data = this->v_ys.GetData();
   double* zs_data = this->v_zs.GetData();

   // set gaps and pressure arrays. Note that for this test 
   // the length of the slave nodes array is the same as the master, 
   // which means that it is the total number of nodes in the whole 
   // mesh
   double* gaps, * pressures;
   int numTotalNodes = this->lengthSlaveNodes;
   gaps = new double[ numTotalNodes ];
   pressures = new double[ numTotalNodes ];

   // initialize arrays
   for (int i=0; i<numTotalNodes; ++i)
   {
      gaps[i] = 0.;
      pressures[i] = 1.;
   }

   // initialize
   int err = Initialize( 3, true );

   // setup simple coupling
   SimpleCouplingSetup( 3,
                        tribol::MORTAR_WEIGHTS,
                        this->numMasterCells,
                        this->lengthMasterNodes,
                        ixm_data,
                        xm_data,
                        ym_data,
                        zm_data,
                        this->numSlaveCells,
                        this->lengthSlaveNodes, 
                        ixs_data,
                        xs_data,
                        ys_data,
                        zs_data,
                        1.e-3,
                        gaps,
                        pressures);

   double dt = 1.0;
   err = Update( dt );

   EXPECT_EQ(err, 0);

   tribol::CouplingSchemeManager& couplingSchemeManager = 
         tribol::CouplingSchemeManager::getInstance();
  
   tribol::CouplingScheme* couplingScheme = couplingSchemeManager.getCoupling( 0 );
   TestMortarWeights( couplingScheme, 2.256, 1.e-3 );

}

TEST_F( MortarSparseWtsTest, mortar_sphere_offset )
{

   // read data sets for mesh
   std::string ixm_file = axom_fs::joinPath(TRIBOL_DATA_DIR,"mortar_sphere_offset/ixm.txt");
   std::string ixs_file = axom_fs::joinPath(TRIBOL_DATA_DIR,"mortar_sphere_offset/ixs.txt");
   std::string xm_file =  axom_fs::joinPath(TRIBOL_DATA_DIR,"mortar_sphere_offset/xm.txt");
   std::string ym_file =  axom_fs::joinPath(TRIBOL_DATA_DIR,"mortar_sphere_offset/ym.txt");
   std::string zm_file =  axom_fs::joinPath(TRIBOL_DATA_DIR,"mortar_sphere_offset/zm.txt");
   std::string xs_file =  axom_fs::joinPath(TRIBOL_DATA_DIR,"mortar_sphere_offset/xs.txt");
   std::string ys_file =  axom_fs::joinPath(TRIBOL_DATA_DIR,"mortar_sphere_offset/ys.txt");
   std::string zs_file =  axom_fs::joinPath(TRIBOL_DATA_DIR,"mortar_sphere_offset/zs.txt");

   this->lengthMasterConn  = 4*507;
   this->lengthSlaveConn   = 4*48;
   this->lengthMasterNodes = 2918;
   this->lengthSlaveNodes  = 2918;
   this->numMasterCells = 507;
   this->numSlaveCells = 48;

   std::ifstream i_ixm( ixm_file ); 
   std::ifstream i_ixs( ixs_file );
   std::ifstream i_xm( xm_file );
   std::ifstream i_ym( ym_file );
   std::ifstream i_zm( zm_file );
   std::ifstream i_xs( xs_file );
   std::ifstream i_ys( ys_file );
   std::ifstream i_zs( zs_file );

   this->v_ixm.SetSize( this->lengthMasterConn );
   this->v_ixs.SetSize( this->lengthSlaveConn );

   this->v_ixm.Load( i_ixm, 1 );
   this->v_ixs.Load( i_ixs, 1 );
   this->v_xm.Load( i_xm, this->lengthMasterNodes );
   this->v_ym.Load( i_ym, this->lengthMasterNodes );
   this->v_zm.Load( i_zm, this->lengthMasterNodes );
   this->v_xs.Load( i_xs, this->lengthSlaveNodes );
   this->v_ys.Load( i_ys, this->lengthSlaveNodes );
   this->v_zs.Load( i_zs, this->lengthSlaveNodes );

   i_ixm.close();
   i_ixs.close();
   i_xm.close();
   i_ym.close();
   i_zm.close();
   i_xs.close();
   i_ys.close();
   i_zs.close();

   SLIC_INFO("After loading mesh data and constructing mfem vectors.");

   // get pointers to mfem vector data
   int* ixm_data   = this->v_ixm.GetData();
   int* ixs_data   = this->v_ixs.GetData();
   double* xm_data = this->v_xm.GetData();
   double* ym_data = this->v_ym.GetData();
   double* zm_data = this->v_zm.GetData();
   double* xs_data = this->v_xs.GetData();
   double* ys_data = this->v_ys.GetData();
   double* zs_data = this->v_zs.GetData();

   // set gaps and pressure arrays. Note that for this test 
   // the length of the slave nodes array is the same as the master, 
   // which means that it is the total number of nodes in the whole 
   // mesh
   double* gaps, * pressures;
   int numTotalNodes = this->lengthSlaveNodes;
   gaps = new double[ numTotalNodes ];
   pressures = new double[ numTotalNodes ];

   // initialize arrays
   for (int i=0; i<numTotalNodes; ++i)
   {
      gaps[i] = 0.;
      pressures[i] = 1.;
   }

   // initialize
   int err = Initialize( 3, true );

   // setup simple coupling
   SimpleCouplingSetup( 3,
                        tribol::MORTAR_WEIGHTS,
                        this->numMasterCells,
                        this->lengthMasterNodes,
                        ixm_data,
                        xm_data,
                        ym_data,
                        zm_data,
                        this->numSlaveCells,
                        this->lengthSlaveNodes, 
                        ixs_data,
                        xs_data,
                        ys_data,
                        zs_data,
                        1.e-3,
                        gaps,
                        pressures);

   double dt = 1.0;
   err = Update( dt );

   EXPECT_EQ(err, 0);

   tribol::CouplingSchemeManager& couplingSchemeManager = 
         tribol::CouplingSchemeManager::getInstance();
  
   tribol::CouplingScheme* couplingScheme = couplingSchemeManager.getCoupling( 0 );
   TestMortarWeights( couplingScheme, 2.260, 1.e-1 );

}

TEST_F( MortarSparseWtsTest, mortar_one_seg_rotated )
{

   // read data sets for mesh
   std::string ixm_file = axom_fs::joinPath(TRIBOL_DATA_DIR,"mortar_one_seg_rotated_square/ixm.txt");
   std::string ixs_file = axom_fs::joinPath(TRIBOL_DATA_DIR,"mortar_one_seg_rotated_square/ixs.txt");
   std::string xm_file =  axom_fs::joinPath(TRIBOL_DATA_DIR,"mortar_one_seg_rotated_square/xm.txt");
   std::string ym_file =  axom_fs::joinPath(TRIBOL_DATA_DIR,"mortar_one_seg_rotated_square/ym.txt");
   std::string zm_file =  axom_fs::joinPath(TRIBOL_DATA_DIR,"mortar_one_seg_rotated_square/zm.txt");
   std::string xs_file =  axom_fs::joinPath(TRIBOL_DATA_DIR,"mortar_one_seg_rotated_square/xs.txt");
   std::string ys_file =  axom_fs::joinPath(TRIBOL_DATA_DIR,"mortar_one_seg_rotated_square/ys.txt");
   std::string zs_file =  axom_fs::joinPath(TRIBOL_DATA_DIR,"mortar_one_seg_rotated_square/zs.txt");

   this->lengthMasterConn  = 4;
   this->lengthSlaveConn   = 4;
   this->lengthMasterNodes = 16;
   this->lengthSlaveNodes  = 16;
   this->numMasterCells = 1;
   this->numSlaveCells = 1;

   std::ifstream i_ixm( ixm_file ); 
   std::ifstream i_ixs( ixs_file );
   std::ifstream i_xm( xm_file );
   std::ifstream i_ym( ym_file );
   std::ifstream i_zm( zm_file );
   std::ifstream i_xs( xs_file );
   std::ifstream i_ys( ys_file );
   std::ifstream i_zs( zs_file );

   this->v_ixm.SetSize( this->lengthMasterConn );
   this->v_ixs.SetSize( this->lengthSlaveConn );

   this->v_ixm.Load( i_ixm, 1 );
   this->v_ixs.Load( i_ixs, 1 );
   this->v_xm.Load( i_xm, this->lengthMasterNodes );
   this->v_ym.Load( i_ym, this->lengthMasterNodes );
   this->v_zm.Load( i_zm, this->lengthMasterNodes );
   this->v_xs.Load( i_xs, this->lengthSlaveNodes );
   this->v_ys.Load( i_ys, this->lengthSlaveNodes );
   this->v_zs.Load( i_zs, this->lengthSlaveNodes );

   i_ixm.close();
   i_ixs.close();
   i_xm.close();
   i_ym.close();
   i_zm.close();
   i_xs.close();
   i_ys.close();
   i_zs.close();

   SLIC_INFO("After loading mesh data and constructing mfem vectors.");

   // get pointers to mfem vector data
   int* ixm_data   = this->v_ixm.GetData();
   int* ixs_data   = this->v_ixs.GetData();
   double* xm_data = this->v_xm.GetData();
   double* ym_data = this->v_ym.GetData();
   double* zm_data = this->v_zm.GetData();
   double* xs_data = this->v_xs.GetData();
   double* ys_data = this->v_ys.GetData();
   double* zs_data = this->v_zs.GetData();

   // set gaps and pressure arrays. Note that for this test 
   // the length of the slave nodes array is the same as the master, 
   // which means that it is the total number of nodes in the whole 
   // mesh
   double* gaps, * pressures;
   int numTotalNodes = this->lengthSlaveNodes;
   gaps = new double[ numTotalNodes ];
   pressures = new double[ numTotalNodes ];

   // initialize arrays
   for (int i=0; i<numTotalNodes; ++i)
   {
      gaps[i] = 0.;
      pressures[i] = 1.;
   }

   // initialize
   int err = Initialize( 3, true );

   // setup simple coupling
   SimpleCouplingSetup( 3,
                        tribol::MORTAR_WEIGHTS,
                        this->numMasterCells,
                        this->lengthMasterNodes,
                        ixm_data,
                        xm_data,
                        ym_data,
                        zm_data,
                        this->numSlaveCells,
                        this->lengthSlaveNodes, 
                        ixs_data,
                        xs_data,
                        ys_data,
                        zs_data,
                        1.e-3,
                        gaps,
                        pressures);

   double dt = 1.0;
   err = Update( dt );

   EXPECT_EQ(err, 0);

   tribol::CouplingSchemeManager& couplingSchemeManager = 
         tribol::CouplingSchemeManager::getInstance();
  
   tribol::CouplingScheme* couplingScheme = couplingSchemeManager.getCoupling( 0 );
   TestMortarWeights( couplingScheme, 20., 1.e-3 );

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
