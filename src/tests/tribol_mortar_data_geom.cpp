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

/*!
 * Test fixture class with some setup necessary to test
 * the _MORTAR_ computational geometry implementation in Tribol
 */
class MortarGeomTest : public ::testing::Test
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

TEST_F( MortarGeomTest, mortar_good_patch )
{

   // read data sets for mesh
   std::string ixm_file = axom_fs::joinPath(TRIBOL_DATA_DIR,"mortar_good_patch/ixm.txt"); 
   std::string ixs_file = axom_fs::joinPath(TRIBOL_DATA_DIR,"mortar_good_patch/ixs.txt");
   std::string xm_file =  axom_fs::joinPath(TRIBOL_DATA_DIR,"mortar_good_patch/xm.txt");
   std::string ym_file =  axom_fs::joinPath(TRIBOL_DATA_DIR,"mortar_good_patch/ym.txt");
   std::string zm_file =  axom_fs::joinPath(TRIBOL_DATA_DIR,"mortar_good_patch/zm.txt");
   std::string xs_file =  axom_fs::joinPath(TRIBOL_DATA_DIR,"mortar_good_patch/xs.txt");
   std::string ys_file =  axom_fs::joinPath(TRIBOL_DATA_DIR,"mortar_good_patch/ys.txt");
   std::string zs_file =  axom_fs::joinPath(TRIBOL_DATA_DIR,"mortar_good_patch/zs.txt");

   this->lengthMasterConn  = 64;
   this->lengthSlaveConn   = 36;
   this->lengthMasterNodes = 189;
   this->lengthSlaveNodes  = 189;
   this->numMasterCells = 16;
   this->numSlaveCells = 9;

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

   EXPECT_EQ( couplingScheme->getNumActivePairs(), 36 );
}

TEST_F( MortarGeomTest, mortar_bad_patch )
{

   // read data sets for mesh
   std::string ixm_file = axom_fs::joinPath(TRIBOL_DATA_DIR,"mortar_bad_patch/ixm.txt"); 
   std::string ixs_file = axom_fs::joinPath(TRIBOL_DATA_DIR,"mortar_bad_patch/ixs.txt");
   std::string xm_file =  axom_fs::joinPath(TRIBOL_DATA_DIR,"mortar_bad_patch/xm.txt");
   std::string ym_file =  axom_fs::joinPath(TRIBOL_DATA_DIR,"mortar_bad_patch/ym.txt");
   std::string zm_file =  axom_fs::joinPath(TRIBOL_DATA_DIR,"mortar_bad_patch/zm.txt");
   std::string xs_file =  axom_fs::joinPath(TRIBOL_DATA_DIR,"mortar_bad_patch/xs.txt");
   std::string ys_file =  axom_fs::joinPath(TRIBOL_DATA_DIR,"mortar_bad_patch/ys.txt");
   std::string zs_file =  axom_fs::joinPath(TRIBOL_DATA_DIR,"mortar_bad_patch/zs.txt");

   this->lengthMasterConn  = 64;
   this->lengthSlaveConn   = 36;
   this->lengthMasterNodes = 189;
   this->lengthSlaveNodes  = 189;
   this->numMasterCells = 16;
   this->numSlaveCells = 9;

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

   EXPECT_EQ( couplingScheme->getNumActivePairs(), 36 );

}

TEST_F( MortarGeomTest, mortar_ironing )
{

   // read data sets for mesh
   std::string ixm_file = axom_fs::joinPath(TRIBOL_DATA_DIR,"mortar_ironing_3/ixm.txt");
   std::string ixs_file = axom_fs::joinPath(TRIBOL_DATA_DIR,"mortar_ironing_3/ixs.txt");
   std::string xm_file =  axom_fs::joinPath(TRIBOL_DATA_DIR,"mortar_ironing_3/xm.txt");
   std::string ym_file =  axom_fs::joinPath(TRIBOL_DATA_DIR,"mortar_ironing_3/ym.txt");
   std::string zm_file =  axom_fs::joinPath(TRIBOL_DATA_DIR,"mortar_ironing_3/zm.txt");
   std::string xs_file =  axom_fs::joinPath(TRIBOL_DATA_DIR,"mortar_ironing_3/xs.txt");
   std::string ys_file =  axom_fs::joinPath(TRIBOL_DATA_DIR,"mortar_ironing_3/ys.txt");
   std::string zs_file =  axom_fs::joinPath(TRIBOL_DATA_DIR,"mortar_ironing_3/zs.txt");

   this->lengthMasterConn  = 480;
   this->lengthSlaveConn   = 160;
   this->lengthMasterNodes = 696;
   this->lengthSlaveNodes  = 696;
   this->numMasterCells = 120;
   this->numSlaveCells = 40;

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

   EXPECT_EQ( err, 0 );

   tribol::CouplingSchemeManager& couplingSchemeManager = 
         tribol::CouplingSchemeManager::getInstance();
  
   tribol::CouplingScheme* couplingScheme = couplingSchemeManager.getCoupling( 0 );

   tribol::IndexType const masterId = couplingScheme->getMeshId1();
   //tribol::IndexType const slaveId = couplingScheme->getMeshId2();
   tribol::MeshManager& meshManager = tribol::MeshManager::getInstance();
   tribol::MeshData& masterMesh = meshManager.GetMeshInstance( masterId );

   if (masterMesh.m_sortedSurfaceNodeIds == nullptr)
   {
      masterMesh.sortSurfaceNodeIds();
   }

   //int nodeOffset = masterMesh.m_sortedSurfaceNodeIds[ masterMesh.m_numSurfaceNodes-1 ] + 1;

   int *I = nullptr;
   int *J = nullptr;
   real *wts = nullptr;
   int nOffsets = 0;
   int nNonZeros = 0;
   int csr_err = GetSimpleCouplingCSR( &I, &J, &wts, &nOffsets, &nNonZeros );

   EXPECT_EQ( csr_err, 0 );

   numTotalNodes = static_cast<tribol::MortarData*>( couplingScheme->getMethodData() )->m_numTotalNodes;
   int num_total_active_nodes = 0;
   for (int a=0; a<numTotalNodes; ++a)
   {
      // loop over range of nonzero column entries
      int num_active_cnt = 0;
      for (int b=I[a]; b<I[a+1]; ++b)
      {
         if (num_active_cnt == 0)
         {
            ++num_active_cnt;
            ++num_total_active_nodes;
         }
      } // end loop over nonzero columns, I[a]
   } // end loop over matrix rows

   SLIC_INFO("Total number of ACTIVE slave nodes: " << num_total_active_nodes );

   EXPECT_EQ( num_total_active_nodes, 54 );

}

TEST_F( MortarGeomTest, mortar_ironing_block_sub_mesh )
{
   // this mesh tests if the update occurred without errors. 
   // Specifically, the Tribol issue was that the inverse 
   // isoparametric mapping for this mesh configuration was providing 
   // a (xi,eta) point outside the parent quadrilateral. A check 
   // was added to the InvIso() routine to catch this and error out.

   // read data sets for mesh
   std::string ixm_file = axom_fs::joinPath(TRIBOL_DATA_DIR,"mortar_ironing_sub_mesh/ixm.txt");
   std::string ixs_file = axom_fs::joinPath(TRIBOL_DATA_DIR,"mortar_ironing_sub_mesh/ixs.txt");
   std::string xm_file =  axom_fs::joinPath(TRIBOL_DATA_DIR,"mortar_ironing_sub_mesh/xm.txt");
   std::string ym_file =  axom_fs::joinPath(TRIBOL_DATA_DIR,"mortar_ironing_sub_mesh/ym.txt");
   std::string zm_file =  axom_fs::joinPath(TRIBOL_DATA_DIR,"mortar_ironing_sub_mesh/zm.txt");
   std::string xs_file =  axom_fs::joinPath(TRIBOL_DATA_DIR,"mortar_ironing_sub_mesh/xs.txt");
   std::string ys_file =  axom_fs::joinPath(TRIBOL_DATA_DIR,"mortar_ironing_sub_mesh/ys.txt");
   std::string zs_file =  axom_fs::joinPath(TRIBOL_DATA_DIR,"mortar_ironing_sub_mesh/zs.txt");

   this->lengthMasterConn  = 28;
   this->lengthSlaveConn   = 8;
   this->lengthMasterNodes = 663;
   this->lengthSlaveNodes  = 663;
   this->numMasterCells = 7;
   this->numSlaveCells = 2;

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

   EXPECT_EQ( err, 0 );

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
