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
#include "axom/core.hpp"
#include "axom/slic.hpp"

#ifdef TRIBOL_USE_UMPIRE
// Umpire includes
#include "umpire/ResourceManager.hpp"
#endif

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
 * the computational geometry. This test does not have a specific 
 * check that will make it pass or fail (yet), instead this is 
 * simply used to drive the computational geometry engine 
 * and to interrogate SLIC output printed to screen
 */
class CommonPlaneCompGeomTest : public ::testing::Test
{
   
public:

   mfem::Vector v_x1;
   mfem::Vector v_y1;
   mfem::Vector v_z1;
   mfem::Vector v_x2;
   mfem::Vector v_y2;
   mfem::Vector v_z2;

   mfem::Array<int> conn1;
   mfem::Array<int> conn2;

   mfem::Vector v_fx1;
   mfem::Vector v_fy1;
   mfem::Vector v_fz1;
   mfem::Vector v_fx2;
   mfem::Vector v_fy2;
   mfem::Vector v_fz2;

   int lengthConn1;
   int lengthConn2;
   int lengthNodes1;
   int lengthNodes2;
   int numCells1;
   int numCells2;
   int numNodesPerCell1;
   int numNodesPerCell2;

   bool isMesh1DataSizeSet {false};
   bool isMesh2DataSizeSet {false};

   void SetMesh1DataSize( const int numTotalNodes, const int numNodesPerCell, const int numFaces )
   {
      this->lengthConn1 = numNodesPerCell * numFaces;
      this->lengthNodes1 = numTotalNodes;
      this->numCells1 = numFaces;
      this->numNodesPerCell1 = numNodesPerCell;
      this->isMesh1DataSizeSet = true;
   }

   void SetMesh2DataSize( const int numTotalNodes, const int numNodesPerCell, const int numFaces )
   {
      this->lengthConn2 = numNodesPerCell * numFaces;
      this->lengthNodes2 = numTotalNodes;
      this->numCells2 = numFaces;
      this->numNodesPerCell2 = numNodesPerCell;
      this->isMesh2DataSizeSet = true;
   }

   void LoadFace1VertexDataFiles3D( std::string &face1x_file, std::string &face1y_file, std::string &face1z_file )
   {
      SLIC_ERROR_IF(!isMesh1DataSizeSet, "LoadFace1VertexDataFiles3D(): Test must call SetMesh1DataSize() first.");

      std::ifstream face1x( face1x_file );
      std::ifstream face1y( face1y_file );
      std::ifstream face1z( face1z_file );
 
      this->v_x1.Load( face1x, this->numCells1*this->numNodesPerCell1 );
      this->v_y1.Load( face1y, this->numCells1*this->numNodesPerCell1 );
      this->v_z1.Load( face1z, this->numCells1*this->numNodesPerCell1 );

      face1x.close();
      face1y.close();
      face1z.close();
   }

   void LoadFace2VertexDataFiles3D( std::string &face2x_file, std::string &face2y_file, std::string &face2z_file )
   {
      SLIC_ERROR_IF(!isMesh2DataSizeSet, "LoadFace2VertexDataFiles3D(): Test must call SetMesh2DataSize() first.");

      std::ifstream face2x( face2x_file );
      std::ifstream face2y( face2y_file );
      std::ifstream face2z( face2z_file );

      this->v_x2.Load( face2x, this->numCells2*this->numNodesPerCell2 );
      this->v_y2.Load( face2y, this->numCells2*this->numNodesPerCell2 );
      this->v_z2.Load( face2z, this->numCells2*this->numNodesPerCell2 );

      face2x.close();
      face2y.close();
      face2z.close();
   } 

   // set the size and initialize the force vectors
   void SetAndInitializeForceVectors()
   {
      SLIC_ERROR_IF(!isMesh1DataSizeSet || !isMesh2DataSizeSet, 
                    "SetAndInitializeForceVectors: Test must call SetMesh1DataSize() and SetMesh2DataSize() first.");
      this->v_fx1.SetSize( this->lengthNodes1 );
      this->v_fy1.SetSize( this->lengthNodes1 );
      this->v_fz1.SetSize( this->lengthNodes1 );
      this->v_fx2.SetSize( this->lengthNodes2 );
      this->v_fy2.SetSize( this->lengthNodes2 );
      this->v_fz2.SetSize( this->lengthNodes2 );

      this->v_fx1 = 0.;
      this->v_fy1 = 0.;
      this->v_fz1 = 0.;
      this->v_fx2 = 0.;
      this->v_fy2 = 0.;
      this->v_fz2 = 0.;
   }

   // set and assign the connectivity of mesh 1 from connectivity pointer
   void SetAndAssignConn1( int* const conn1 )
   {
      SLIC_ERROR_IF(!isMesh1DataSizeSet, "SetAndAssignConn1(): Test must call SetMesh1DataSize() first.");
      this->conn1.SetSize( this->numNodesPerCell1 * this->numCells1 );
      this->conn1.Assign(conn1);
   }

   // set and assign the connectivity of mesh 2 from connectivity pointer
   void SetAndAssignConn2( int* const conn2 )
   {
      SLIC_ERROR_IF(!isMesh2DataSizeSet, "SetAndAssignConn2(): Test must call SetMesh2DataSize() first.");
      this->conn2.SetSize( this->numNodesPerCell2 * this->numCells2 );
      this->conn2.Assign(conn2);
   }


protected:

   void SetUp() override
   {
   }

   void TearDown() override
   {
   }

protected:

};

TEST_F( CommonPlaneCompGeomTest, overlap_with_duplicate_vertex_misordering )
{
   // This tests the computational geometry on a single problematic face-pair 
   // that arose in host-code testing
   int dim = 3;
   int numVertsPerFace = 4;
   int numFacesPerMesh = 1;
   int numTotalNodesPerMesh = numVertsPerFace;

   // read data sets for mesh
   std::string face1x_f = axom_fs::joinPath(TRIBOL_DATA_DIR,"common_plane_parallel_1_overlap/face1x_bug.txt");
   std::string face1y_f = axom_fs::joinPath(TRIBOL_DATA_DIR,"common_plane_parallel_1_overlap/face1y_bug.txt");
   std::string face1z_f = axom_fs::joinPath(TRIBOL_DATA_DIR,"common_plane_parallel_1_overlap/face1z_bug.txt");
   std::string face2x_f = axom_fs::joinPath(TRIBOL_DATA_DIR,"common_plane_parallel_1_overlap/face2x_bug.txt");
   std::string face2y_f = axom_fs::joinPath(TRIBOL_DATA_DIR,"common_plane_parallel_1_overlap/face2y_bug.txt");
   std::string face2z_f = axom_fs::joinPath(TRIBOL_DATA_DIR,"common_plane_parallel_1_overlap/face2z_bug.txt");

   SetMesh1DataSize( numTotalNodesPerMesh, numVertsPerFace, numFacesPerMesh );
   SetMesh2DataSize( numTotalNodesPerMesh, numVertsPerFace, numFacesPerMesh );

   LoadFace1VertexDataFiles3D( face1x_f, face1y_f, face1z_f );
   LoadFace2VertexDataFiles3D( face2x_f, face2y_f, face2z_f );

   SetAndInitializeForceVectors();

   // manually setup the connectivity for this single face per side problem
   int face_conn1[numVertsPerFace] = {0,1,2,3};
   int face_conn2[numVertsPerFace] = {0,1,2,3};

   SetAndAssignConn1( &face_conn1[0] );
   SetAndAssignConn2( &face_conn2[0] );

   tribol::CommType problem_comm = TRIBOL_COMM_WORLD;
   tribol::initialize( dim, problem_comm );

   tribol::registerMesh( 0, numFacesPerMesh, numFacesPerMesh*numVertsPerFace, this->conn1.GetData(), (int)(tribol::LINEAR_QUAD),
                         this->v_x1.GetData(), this->v_y1.GetData(), this->v_z1.GetData() );
   tribol::registerMesh( 1, numFacesPerMesh, numFacesPerMesh*numVertsPerFace, this->conn2.GetData(), (int)(tribol::LINEAR_QUAD),
                         this->v_x2.GetData(), this->v_y2.GetData(), this->v_z2.GetData() );

   tribol::registerNodalResponse( 0, this->v_fx1.GetData(), this->v_fy1.GetData(), this->v_fz1.GetData() );
   tribol::registerNodalResponse( 0, this->v_fx2.GetData(), this->v_fy2.GetData(), this->v_fz2.GetData() );

   tribol::setKinematicConstantPenalty( 0, 1. );
   tribol::setKinematicConstantPenalty( 1, 1. );

   tribol::registerCouplingScheme( 0, 0, 1,
                                   tribol::SURFACE_TO_SURFACE,
                                   tribol::AUTO,
                                   tribol::COMMON_PLANE,
                                   tribol::FRICTIONLESS,
                                   tribol::PENALTY,
                                   tribol::BINNING_CARTESIAN_PRODUCT );

   tribol::setPenaltyOptions( 0, tribol::KINEMATIC, tribol::KINEMATIC_CONSTANT );
   tribol::setContactPenFrac(0.5);
   tribol::setContactAreaFrac(1.e-4);

   double dt = 1.;
   int update_err = tribol::update( 1, 1., dt );

   EXPECT_EQ( update_err, 0 );
}


int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

#ifdef TRIBOL_USE_UMPIRE
  umpire::ResourceManager::getInstance();         // initialize umpire's ResouceManager
#endif

  axom::slic::SimpleLogger logger;                // create & initialize logger,
  tribol::SimpleMPIWrapper wrapper(argc, argv);   // initialize and finalize MPI, when applicable

  result = RUN_ALL_TESTS();

  return result;
}
