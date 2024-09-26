// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Tribol Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (MIT)

// Tribol includes
#include "tribol/interface/tribol.hpp"
#include "tribol/utils/TestUtils.hpp"
#include "tribol/utils/Math.hpp"
#include "tribol/common/Parameters.hpp"
#include "tribol/mesh/MethodCouplingData.hpp"
#include "tribol/mesh/CouplingScheme.hpp"
#include "tribol/mesh/InterfacePairs.hpp"
#include "tribol/mesh/MeshData.hpp"
#include "tribol/physics/CommonPlane.hpp"
#include "tribol/geom/GeomUtilities.hpp"

#ifdef TRIBOL_USE_UMPIRE
// Umpire includes
#include "umpire/ResourceManager.hpp"
#endif

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

using RealT = tribol::RealT;

void compareGaps( tribol::CouplingScheme const * cs, 
                  RealT gap, const RealT tol,
                  const char *gapType )
{
   tribol::IndexT const numPairs = cs->getNumActivePairs();
   const auto cs_view = const_cast<tribol::CouplingScheme*>(cs)->getView();

   for (tribol::IndexT cpID = 0; cpID < numPairs; ++cpID)
   {
      auto& plane = cs->getContactPlane(cpID);

      RealT my_gap = 0.;
      if ( std::strcmp( gapType, "kinematic_penetration" ) == 0 ||
           std::strcmp( gapType, "kinematic_separation"  ) == 0 )
      {
         my_gap = plane.m_gap;
      }
      else
      {
         my_gap = plane.m_velGap;
      }

      // gap tolerance per common-plane based calculation (see CouplingScheme::ViewerBase::getCommonPlaneGapTol)
      RealT gap_tol = cs_view.getCommonPlaneGapTol( plane.getCpElementId1(), plane.getCpElementId2() );
                      
      // check gap sense.
      if ( std::strcmp( gapType, "kinematic_penetration" ) == 0  || 
           std::strcmp( gapType, "rate_penetration" ) == 0 )
      { 
         // check that the gap is less than the separation tolerance (interpenetration)
         EXPECT_LE( my_gap, gap_tol );
      }
      else if ( std::strcmp( gapType, "kinematic_separation" ) == 0 || 
                std::strcmp( gapType, "rate_separation" ) == 0 )
      {
         // check that the gap is greater than the separation tolerance (separation)
         EXPECT_GE( my_gap, gap_tol );
      }
      else
      {
         SLIC_ERROR("compareGaps: invalid gapType. " <<
                    "Acceptable types are 'kinematic_penetration', 'kinematic_separation', " <<
                    "'rate_penetration' or 'rate_separation'." ) ;;
      }

      // check diffs
      RealT diff = std::abs( my_gap - gap );
      EXPECT_LE( diff, tol );
   }
} // end compareGaps()

void checkMeshPenalties( tribol::CouplingScheme const * cs,
                         const RealT penalty, const RealT tol, 
                         const char * penaltyType )
{
   const tribol::IndexT mesh_id1 = cs->getMeshId1();
   const tribol::IndexT mesh_id2 = cs->getMeshId2();

   tribol::MeshManager& meshManager = tribol::MeshManager::getInstance();
   tribol::MeshData& mesh1 = meshManager.at( mesh_id1 );
   tribol::MeshData& mesh2 = meshManager.at( mesh_id2 );

   if ( std::strcmp( penaltyType, "constant" ) == 0 )
   {
      RealT penalty_diff_1 = std::abs( mesh1.getElementData().m_penalty_stiffness - penalty );
      RealT penalty_diff_2 = std::abs( mesh2.getElementData().m_penalty_stiffness - penalty );
      EXPECT_LE( penalty_diff_1, tol );
      EXPECT_LE( penalty_diff_2, tol );
   }
   else if ( std::strcmp( penaltyType, "face" ) == 0 )
   {
      // no-op, the face-based penalty is checked in a call to tribol::update()
   }
   else if ( std::strcmp( penaltyType, "constant_rate" ) == 0 )
   {
      RealT penalty_diff_1 = std::abs( mesh1.getElementData().m_rate_penalty_stiffness - penalty );
      RealT penalty_diff_2 = std::abs( mesh2.getElementData().m_rate_penalty_stiffness - penalty );
      EXPECT_LE( penalty_diff_1, tol );
      EXPECT_LE( penalty_diff_2, tol );
   }
   else if ( std::strcmp( penaltyType, "percent_rate" ) == 0 )
   {
      RealT penalty1 = mesh1.getElementData().m_rate_percent_stiffness * mesh1.getElementData().m_penalty_stiffness; 
      RealT penalty2 = mesh2.getElementData().m_rate_percent_stiffness * mesh2.getElementData().m_penalty_stiffness; 
      RealT penalty_diff_1 = std::abs( penalty1 - penalty );
      RealT penalty_diff_2 = std::abs( penalty2 - penalty );
      EXPECT_LE( penalty_diff_1, tol );
      EXPECT_LE( penalty_diff_2, tol );
   }
   else
   {
      SLIC_ERROR("checkMeshPenalties: invalid penaltyType. " << 
                 "only 'constant', 'face', 'constant_rate', or 'percent_rate' accepted. " );
   }

} // end checkMeshPenalties()

void checkPressures( tribol::CouplingScheme const * cs, 
                     RealT pressure, const RealT tol, const char * pressureType = "kinematic"  )
{
   tribol::IndexT const numPairs = cs->getNumActivePairs();

   for (tribol::IndexT cpID = 0; cpID < numPairs; ++cpID)
   {
      auto& plane = cs->getContactPlane(cpID);

      RealT my_pressure = 0.;
      if ( std::strcmp( pressureType, "rate" ) == 0 )
      {
         my_pressure = plane.m_ratePressure;
      }
      else if ( std::strcmp( pressureType, "kinematic" ) == 0 )
      {
         my_pressure = plane.m_pressure;
      }
      else
      {
         SLIC_ERROR( "checkPressures(): invalid pressure type. Supported types are " << 
                     "'kinematic' or 'rate'." );
      }

      // check diffs
      RealT press_diff = std::abs( my_pressure - pressure );
      EXPECT_LE( press_diff, tol );
   }
} // end checkPressures()

// problem specific routine to check the sense of the force. Note, this 
// routine makes implicit use of the knowledge that the outward facing 
// surface unit normals are in the +/- z-direction for mesh 1 and 
// mesh 2, respectively. This is not a general routine for general 
// mesh configurations.
void checkForceSense( tribol::CouplingScheme const * cs, bool isTied = false )
{
   // TODO: get rid of const cast
   const auto mesh1 = const_cast<tribol::CouplingScheme*>(cs)->getMesh1().getView();
   const auto mesh2 = const_cast<tribol::CouplingScheme*>(cs)->getMesh2().getView();

   for (int i=0; i<2; ++i) // loop over meshes
   { 
      auto& mesh = (i==0) ? mesh1 : mesh2;
  
      // loop over faces and nodes
      for (tribol::IndexT kf = 0; kf < mesh.numberOfElements(); ++kf)
      {
         for (tribol::IndexT a = 0; a<mesh.numberOfNodesPerElement(); ++a)
         {
            int node_id = mesh.getGlobalNodeId(kf, a);
            RealT force_mag = tribol::dotProd( mesh.getResponse()[0][ node_id ],
                                              mesh.getResponse()[1][ node_id ], 
                                              mesh.getResponse()[2][ node_id ],
                                              mesh.getElementNormals()[0][ kf ],
                                              mesh.getElementNormals()[1][ kf ],
                                              mesh.getElementNormals()[2][ kf ] );
            if (!isTied) {
               // <= catches interpenetration AND separation
               EXPECT_LE( force_mag, 0. );
            }
            else {
               // no-op, TIED_NORMAL is a special case where we 
               // support all force 'sense' (i.e. tension AND compression)
            }
         }
      }
   }
} // end checkForceSense()

/*!
 * Test fixture class with some setup necessary to test 
 * the COMMON_PLANE + GAP RATE PENALTY implementation 
 */
class CommonPlaneTest : public ::testing::Test
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

TEST_F( CommonPlaneTest, constant_rate_penetration )
{
   this->m_mesh.mortarMeshId = 0;
   this->m_mesh.nonmortarMeshId = 1;

   int nMortarElems = 1; 
   int nElemsXM = nMortarElems;
   int nElemsYM = nMortarElems;
   int nElemsZM = nMortarElems;

   int nNonmortarElems = 1; 
   int nElemsXS = nNonmortarElems;
   int nElemsYS = nNonmortarElems;
   int nElemsZS = nNonmortarElems;

   // mesh bounding box with 0.1 interpenetration gap
   RealT x_min1 = 0.;
   RealT y_min1 = 0.;
   RealT z_min1 = 0.; 
   RealT x_max1 = 1.;
   RealT y_max1 = 1.;
   RealT z_max1 = 1.05;

   RealT x_min2 = 0.;
   RealT y_min2 = 0.;
   RealT z_min2 = 0.95;
   RealT x_max2 = 1.;
   RealT y_max2 = 1.;
   RealT z_max2 = 2.;

   this->m_mesh.setupContactMeshHex( nElemsXM, nElemsYM, nElemsZM,
                                     x_min1, y_min1, z_min1,
                                     x_max1, y_max1, z_max1,
                                     nElemsXS, nElemsYS, nElemsZS,
                                     x_min2, y_min2, z_min2,
                                     x_max2, y_max2, z_max2,
                                     0., 0. );

   // set velocities
   RealT dt = 1.e-3;
   RealT velX1 = 0.;
   RealT velY1 = 0.;
   RealT velZ1 = 1.;
   RealT velX2 = 0.;
   RealT velY2 = 0.;
   RealT velZ2 = -1.;
   this->m_mesh.allocateAndSetVelocities( m_mesh.mortarMeshId, velX1, velY1, velZ1 );
   this->m_mesh.allocateAndSetVelocities( m_mesh.nonmortarMeshId,  velX2, velY2, velZ2 );
   
   // set kinematic gap and rate (velocity) gap penalty parameters
   tribol::TestControlParameters parameters;
   parameters.penalty_ratio = false;
   parameters.const_penalty = 0.75;
   parameters.percent_rate_penalty = false; // default
   parameters.rate_penalty_ratio = 0.0; // default
   parameters.constant_rate_penalty = true;
   parameters.rate_penalty = 1.0;
   parameters.dt = dt;

   int test_mesh_update_err = 
      this->m_mesh.tribolSetupAndUpdate( tribol::COMMON_PLANE, tribol::PENALTY, 
                                         tribol::FRICTIONLESS, tribol::NO_CASE, false, parameters );

   EXPECT_EQ( test_mesh_update_err, 0 );

   tribol::CouplingSchemeManager& couplingSchemeManager = 
         tribol::CouplingSchemeManager::getInstance();
  
   tribol::CouplingScheme* couplingScheme = 
      &couplingSchemeManager.at( 0 );

   // check mesh rate penalties
   RealT penalty = 0.;
   if (parameters.constant_rate_penalty)
   {
      penalty = parameters.rate_penalty;
      checkMeshPenalties( couplingScheme, penalty, 1.E-8, "constant_rate" );
   }
   else if (parameters.percent_rate_penalty)
   {
      penalty = parameters.rate_penalty_ratio * parameters.const_penalty;
      checkMeshPenalties( couplingScheme, penalty, 1.E-8, "percent_rate" );
   }

   // check the gaps, pressures, and force sense
   //RealT gap = z_min2 - z_max1;
   RealT rate_gap = velZ2 - velZ1;
   RealT pressure = (rate_gap < 0.) ? penalty * rate_gap : 0.;
   compareGaps( couplingScheme, rate_gap, 1.E-8, "rate_penetration" );
   checkPressures( couplingScheme, pressure, 1.E-8, "rate" );
   checkForceSense( couplingScheme ); // note: the kinematic and rate contributions are not separated

   tribol::finalize();

} // end test 'rate_penetration'

TEST_F( CommonPlaneTest, constant_rate_separation )
{
   this->m_mesh.mortarMeshId = 0;
   this->m_mesh.nonmortarMeshId = 1;

   int nMortarElems = 1; 
   int nElemsXM = nMortarElems;
   int nElemsYM = nMortarElems;
   int nElemsZM = nMortarElems;

   int nNonmortarElems = 1; 
   int nElemsXS = nNonmortarElems;
   int nElemsYS = nNonmortarElems;
   int nElemsZS = nNonmortarElems;

   // mesh bounding box with 0.1 interpenetration kinematic gap
   RealT x_min1 = 0.;
   RealT y_min1 = 0.;
   RealT z_min1 = 0.; 
   RealT x_max1 = 1.;
   RealT y_max1 = 1.;
   RealT z_max1 = 1.05;

   RealT x_min2 = 0.;
   RealT y_min2 = 0.;
   RealT z_min2 = 0.95;
   RealT x_max2 = 1.;
   RealT y_max2 = 1.;
   RealT z_max2 = 2.;

   this->m_mesh.setupContactMeshHex( nElemsXM, nElemsYM, nElemsZM,
                                     x_min1, y_min1, z_min1,
                                     x_max1, y_max1, z_max1,
                                     nElemsXS, nElemsYS, nElemsZS,
                                     x_min2, y_min2, z_min2,
                                     x_max2, y_max2, z_max2,
                                     0., 0. );

   // set velocities such that the blocks have equal and opposite 'separation' velocities
   // Note: this should not trigger a gap-rate contribution.
   RealT dt = 1.e-3;
   RealT velX1 = 0.;
   RealT velY1 = 0.;
   RealT velZ1 = -1.;
   RealT velX2 = 0.;
   RealT velY2 = 0.;
   RealT velZ2 = 1.;
   this->m_mesh.allocateAndSetVelocities( m_mesh.mortarMeshId, velX1, velY1, velZ1 );
   this->m_mesh.allocateAndSetVelocities( m_mesh.nonmortarMeshId,  velX2, velY2, velZ2 );
   
   // set kinematic gap and rate (velocity) gap penalty parameters
   tribol::TestControlParameters parameters;
   parameters.penalty_ratio = false;
   parameters.const_penalty = 0.75;
   parameters.percent_rate_penalty = false; // default
   parameters.rate_penalty_ratio = 0.0; // default
   parameters.constant_rate_penalty = true;
   parameters.rate_penalty = 1.0;
   parameters.dt = dt;

   int test_mesh_update_err = 
      this->m_mesh.tribolSetupAndUpdate( tribol::COMMON_PLANE, tribol::PENALTY, 
                                         tribol::FRICTIONLESS, tribol::NO_CASE, false, parameters );

   EXPECT_EQ( test_mesh_update_err, 0 );

   tribol::CouplingSchemeManager& couplingSchemeManager = 
         tribol::CouplingSchemeManager::getInstance();
  
   tribol::CouplingScheme* couplingScheme = 
      &couplingSchemeManager.at( 0 );

   // check mesh rate penalties
   RealT penalty = 0.;
   if (parameters.constant_rate_penalty)
   {
      penalty = parameters.rate_penalty;
      checkMeshPenalties( couplingScheme, penalty, 1.E-8, "constant_rate" );
   }
   else if (parameters.percent_rate_penalty)
   {
      penalty = parameters.rate_penalty_ratio * parameters.const_penalty;
      checkMeshPenalties( couplingScheme, penalty, 1.E-8, "percent_rate" );
   }

   // check the gaps, pressures and force sense
   //RealT gap = z_min2 - z_max1;
   RealT rate_gap = velZ2 - velZ1;
   RealT pressure = (rate_gap < 0.) ? penalty * rate_gap : 0.;
   compareGaps( couplingScheme, rate_gap, 1.E-8, "rate_separation" );
   checkPressures( couplingScheme, pressure, 1.E-8, "rate" );
   checkForceSense( couplingScheme ); // note: the kinematic and rate contributions aren't separated

   tribol::finalize();

} // end test 'rate_separation'

TEST_F( CommonPlaneTest, no_gap_constant_rate_penetration )
{
   this->m_mesh.mortarMeshId = 0;
   this->m_mesh.nonmortarMeshId = 1;

   int nMortarElems = 1; 
   int nElemsXM = nMortarElems;
   int nElemsYM = nMortarElems;
   int nElemsZM = nMortarElems;

   int nNonmortarElems = 1; 
   int nElemsXS = nNonmortarElems;
   int nElemsYS = nNonmortarElems;
   int nElemsZS = nNonmortarElems;

   // mesh bounding box with zero gap.
   RealT x_min1 = 0.;
   RealT y_min1 = 0.;
   RealT z_min1 = 0.; 
   RealT x_max1 = 1.;
   RealT y_max1 = 1.;
   RealT z_max1 = 1.;

   RealT x_min2 = 0.;
   RealT y_min2 = 0.;
   RealT z_min2 = 1.;
   RealT x_max2 = 1.;
   RealT y_max2 = 1.;
   RealT z_max2 = 2.;

   this->m_mesh.setupContactMeshHex( nElemsXM, nElemsYM, nElemsZM,
                                     x_min1, y_min1, z_min1,
                                     x_max1, y_max1, z_max1,
                                     nElemsXS, nElemsYS, nElemsZS,
                                     x_min2, y_min2, z_min2,
                                     x_max2, y_max2, z_max2,
                                     0., 0. );

   // set velocities such that further interpenetration would occur
   RealT dt = 1.e-3;
   RealT velX1 = 0.;
   RealT velY1 = 0.;
   RealT velZ1 = 1.;
   RealT velX2 = 0.;
   RealT velY2 = 0.;
   RealT velZ2 = -1.;
   this->m_mesh.allocateAndSetVelocities( m_mesh.mortarMeshId, velX1, velY1, velZ1 );
   this->m_mesh.allocateAndSetVelocities( m_mesh.nonmortarMeshId,  velX2, velY2, velZ2 );
   
   // set kinematic gap and rate (velocity) gap penalty parameters
   tribol::TestControlParameters parameters;
   parameters.penalty_ratio = false;
   parameters.const_penalty = 0.75;
   parameters.percent_rate_penalty = false; // default
   parameters.rate_penalty_ratio = 0.0; // default
   parameters.constant_rate_penalty = true;
   parameters.rate_penalty = 1.0;
   parameters.dt = dt;

   int test_mesh_update_err = 
      this->m_mesh.tribolSetupAndUpdate( tribol::COMMON_PLANE, tribol::PENALTY, 
                                         tribol::FRICTIONLESS, tribol::NO_CASE, false, parameters );

   EXPECT_EQ( test_mesh_update_err, 0 );

   tribol::CouplingSchemeManager& couplingSchemeManager = 
         tribol::CouplingSchemeManager::getInstance();
  
   tribol::CouplingScheme* couplingScheme = 
      &couplingSchemeManager.at( 0 );

   // check mesh rate penalties
   RealT penalty = 0.;
   if (parameters.constant_rate_penalty)
   {
      penalty = parameters.rate_penalty;
      checkMeshPenalties( couplingScheme, penalty, 1.E-8, "constant_rate" );
   }
   else if (parameters.percent_rate_penalty)
   {
      penalty = parameters.rate_penalty_ratio * parameters.const_penalty;
      checkMeshPenalties( couplingScheme, penalty, 1.E-8, "percent_rate" );
   }

   // check to make sure zero rate pressure
   RealT gap = z_min2 - z_max1;
   RealT rate_gap = velZ2 - velZ1;
   RealT pressure = (gap < 0. && rate_gap < 0.) ? penalty * rate_gap : 0.;
   checkPressures( couplingScheme, pressure, 1.E-8, "rate" );

   tribol::finalize();
}

TEST_F( CommonPlaneTest, percent_rate_penetration )
{
   this->m_mesh.mortarMeshId = 0;
   this->m_mesh.nonmortarMeshId = 1;

   int nMortarElems = 1; 
   int nElemsXM = nMortarElems;
   int nElemsYM = nMortarElems;
   int nElemsZM = nMortarElems;

   int nNonmortarElems = 1; 
   int nElemsXS = nNonmortarElems;
   int nElemsYS = nNonmortarElems;
   int nElemsZS = nNonmortarElems;

   // mesh bounding box with 0.1 interpenetration gap
   RealT x_min1 = 0.;
   RealT y_min1 = 0.;
   RealT z_min1 = 0.; 
   RealT x_max1 = 1.;
   RealT y_max1 = 1.;
   RealT z_max1 = 1.05;

   RealT x_min2 = 0.;
   RealT y_min2 = 0.;
   RealT z_min2 = 0.95;
   RealT x_max2 = 1.;
   RealT y_max2 = 1.;
   RealT z_max2 = 2.;

   this->m_mesh.setupContactMeshHex( nElemsXM, nElemsYM, nElemsZM,
                                     x_min1, y_min1, z_min1,
                                     x_max1, y_max1, z_max1,
                                     nElemsXS, nElemsYS, nElemsZS,
                                     x_min2, y_min2, z_min2,
                                     x_max2, y_max2, z_max2,
                                     0., 0. );

   // set velocities such that further penetration would occur
   RealT dt = 1.e-3;
   RealT velX1 = 0.;
   RealT velY1 = 0.;
   RealT velZ1 = 1.;
   RealT velX2 = 0.;
   RealT velY2 = 0.;
   RealT velZ2 = -1.;
   this->m_mesh.allocateAndSetVelocities( m_mesh.mortarMeshId, velX1, velY1, velZ1 );
   this->m_mesh.allocateAndSetVelocities( m_mesh.nonmortarMeshId,  velX2, velY2, velZ2 );
   
   // set kinematic gap and rate (velocity) gap penalty parameters
   tribol::TestControlParameters parameters;
   parameters.penalty_ratio = false;
   parameters.const_penalty = 0.75;
   parameters.percent_rate_penalty = true;
   parameters.rate_penalty_ratio = 0.20;
   parameters.constant_rate_penalty = false;
   parameters.rate_penalty = 1.0;
   parameters.dt = dt;

   int test_mesh_update_err = 
      this->m_mesh.tribolSetupAndUpdate( tribol::COMMON_PLANE, tribol::PENALTY, 
                                         tribol::FRICTIONLESS, tribol::NO_CASE, false, parameters );

   EXPECT_EQ( test_mesh_update_err, 0 );

   tribol::CouplingSchemeManager& couplingSchemeManager = 
         tribol::CouplingSchemeManager::getInstance();
  
   tribol::CouplingScheme* couplingScheme = 
      &couplingSchemeManager.at( 0 );

   // check mesh rate penalties
   RealT penalty = 0.;
   if (parameters.constant_rate_penalty)
   {
      penalty = parameters.rate_penalty;
      checkMeshPenalties( couplingScheme, penalty, 1.E-8, "constant_rate" );
   }
   else if (parameters.percent_rate_penalty)
   {
      penalty = parameters.rate_penalty_ratio * parameters.const_penalty;
      checkMeshPenalties( couplingScheme, penalty, 1.E-8, "percent_rate" );
   }

   // check the gaps, pressures, and force sense
   //RealT gap = z_min2 - z_max1;
   RealT rate_gap = velZ2 - velZ1;
   RealT stiffness = tribol::ComputePenaltyStiffnessPerArea( penalty, penalty );
   RealT pressure = (rate_gap < 0.) ? stiffness * rate_gap : 0.;
   compareGaps( couplingScheme, rate_gap, 1.E-8, "rate_penetration" );
   checkPressures( couplingScheme, pressure, 1.E-8, "rate" );
   checkForceSense( couplingScheme ); // note: the kinematic and rate contributions are not separated

   tribol::finalize();

} // end test 'rate_penetration'

TEST_F( CommonPlaneTest, percent_rate_separation )
{
   this->m_mesh.mortarMeshId = 0;
   this->m_mesh.nonmortarMeshId = 1;

   int nMortarElems = 1; 
   int nElemsXM = nMortarElems;
   int nElemsYM = nMortarElems;
   int nElemsZM = nMortarElems;

   int nNonmortarElems = 1; 
   int nElemsXS = nNonmortarElems;
   int nElemsYS = nNonmortarElems;
   int nElemsZS = nNonmortarElems;

   // mesh bounding box with 0.1 interpenetration gap
   RealT x_min1 = 0.;
   RealT y_min1 = 0.;
   RealT z_min1 = 0.; 
   RealT x_max1 = 1.;
   RealT y_max1 = 1.;
   RealT z_max1 = 1.05;

   RealT x_min2 = 0.;
   RealT y_min2 = 0.;
   RealT z_min2 = 0.95;
   RealT x_max2 = 1.;
   RealT y_max2 = 1.;
   RealT z_max2 = 2.;

   this->m_mesh.setupContactMeshHex( nElemsXM, nElemsYM, nElemsZM,
                                     x_min1, y_min1, z_min1,
                                     x_max1, y_max1, z_max1,
                                     nElemsXS, nElemsYS, nElemsZS,
                                     x_min2, y_min2, z_min2,
                                     x_max2, y_max2, z_max2,
                                     0., 0. );

   // set velocities such that equal and opposite 'separation' velocities would 
   // lead to eventual separation of the blocks. This should not trigger a gap-rate
   // contribution.
   RealT dt = 1.e-3;
   RealT velX1 = 0.;
   RealT velY1 = 0.;
   RealT velZ1 = -1.;
   RealT velX2 = 0.;
   RealT velY2 = 0.;
   RealT velZ2 = 1.;
   this->m_mesh.allocateAndSetVelocities( m_mesh.mortarMeshId, velX1, velY1, velZ1 );
   this->m_mesh.allocateAndSetVelocities( m_mesh.nonmortarMeshId,  velX2, velY2, velZ2 );
   
   // set kinematic gap and rate (velocity) gap penalty parameters
   tribol::TestControlParameters parameters;
   parameters.penalty_ratio = false;
   parameters.const_penalty = 0.75;
   parameters.percent_rate_penalty = true; 
   parameters.rate_penalty_ratio = 0.20;
   parameters.constant_rate_penalty = false;
   parameters.rate_penalty = 1.0;
   parameters.dt = dt;

   int test_mesh_update_err = 
      this->m_mesh.tribolSetupAndUpdate( tribol::COMMON_PLANE, tribol::PENALTY, 
                                         tribol::FRICTIONLESS, tribol::NO_CASE, false, parameters );

   EXPECT_EQ( test_mesh_update_err, 0 );

   tribol::CouplingSchemeManager& couplingSchemeManager = 
         tribol::CouplingSchemeManager::getInstance();
  
   tribol::CouplingScheme* couplingScheme = 
      &couplingSchemeManager.at( 0 );

   // check mesh rate penalties
   RealT penalty = 0.;
   if (parameters.constant_rate_penalty)
   {
      penalty = parameters.rate_penalty;
      checkMeshPenalties( couplingScheme, penalty, 1.E-8, "constant_rate" );
   }
   else if (parameters.percent_rate_penalty)
   {
      penalty = parameters.rate_penalty_ratio * parameters.const_penalty;
      checkMeshPenalties( couplingScheme, penalty, 1.E-8, "percent_rate" );
   }

   // check the gaps, pressures and force sense
   //RealT gap = z_min2 - z_max1;
   RealT rate_gap = velZ2 - velZ1;
   RealT pressure = (rate_gap < 0.) ? penalty * rate_gap : 0.;
   compareGaps( couplingScheme, rate_gap, 1.E-8, "rate_separation" );
   checkPressures( couplingScheme, pressure, 1.E-8, "rate" );
   checkForceSense( couplingScheme ); // note: the kinematic and rate contributions aren't separated

   tribol::finalize();

} // end test 'rate_separation'

TEST_F( CommonPlaneTest, no_gap_percent_rate_penetration )
{
   this->m_mesh.mortarMeshId = 0;
   this->m_mesh.nonmortarMeshId = 1;

   int nMortarElems = 1; 
   int nElemsXM = nMortarElems;
   int nElemsYM = nMortarElems;
   int nElemsZM = nMortarElems;

   int nNonmortarElems = 1; 
   int nElemsXS = nNonmortarElems;
   int nElemsYS = nNonmortarElems;
   int nElemsZS = nNonmortarElems;

   // mesh bounding box with 0.1 interpenetration gap
   RealT x_min1 = 0.;
   RealT y_min1 = 0.;
   RealT z_min1 = 0.; 
   RealT x_max1 = 1.;
   RealT y_max1 = 1.;
   RealT z_max1 = 1.;

   RealT x_min2 = 0.;
   RealT y_min2 = 0.;
   RealT z_min2 = 1.;
   RealT x_max2 = 1.;
   RealT y_max2 = 1.;
   RealT z_max2 = 2.;

   this->m_mesh.setupContactMeshHex( nElemsXM, nElemsYM, nElemsZM,
                                     x_min1, y_min1, z_min1,
                                     x_max1, y_max1, z_max1,
                                     nElemsXS, nElemsYS, nElemsZS,
                                     x_min2, y_min2, z_min2,
                                     x_max2, y_max2, z_max2,
                                     0., 0. );

   // set velocities such that further interpenetration would occur
   RealT dt = 1.e-3;
   RealT velX1 = 0.;
   RealT velY1 = 0.;
   RealT velZ1 = 1.;
   RealT velX2 = 0.;
   RealT velY2 = 0.;
   RealT velZ2 = -1.;
   this->m_mesh.allocateAndSetVelocities( m_mesh.mortarMeshId, velX1, velY1, velZ1 );
   this->m_mesh.allocateAndSetVelocities( m_mesh.nonmortarMeshId,  velX2, velY2, velZ2 );
   
   // set kinematic gap and rate (velocity) gap penalty parameters
   tribol::TestControlParameters parameters;
   parameters.penalty_ratio = false;
   parameters.const_penalty = 0.75;
   parameters.percent_rate_penalty = true; 
   parameters.rate_penalty_ratio = 0.2;
   parameters.constant_rate_penalty = false;
   parameters.rate_penalty = 1.0;
   parameters.dt = dt;

   int test_mesh_update_err = 
      this->m_mesh.tribolSetupAndUpdate( tribol::COMMON_PLANE, tribol::PENALTY, 
                                         tribol::FRICTIONLESS, tribol::NO_CASE, false, parameters );

   EXPECT_EQ( test_mesh_update_err, 0 );

   tribol::CouplingSchemeManager& couplingSchemeManager = 
         tribol::CouplingSchemeManager::getInstance();
  
   tribol::CouplingScheme* couplingScheme = 
      &couplingSchemeManager.at( 0 );

   // check mesh rate penalties
   RealT penalty = 0.;
   if (parameters.constant_rate_penalty)
   {
      penalty = parameters.rate_penalty;
      checkMeshPenalties( couplingScheme, penalty, 1.E-8, "constant_rate" );
   }
   else if (parameters.percent_rate_penalty)
   {
      penalty = parameters.rate_penalty_ratio * parameters.const_penalty;
      checkMeshPenalties( couplingScheme, penalty, 1.E-8, "percent_rate" );
   }

   // check to make sure zero rate pressure
   RealT gap = z_min2 - z_max1;
   RealT rate_gap = velZ2 - velZ1;
   RealT pressure = (gap < 0. && rate_gap < 0.) ? penalty * rate_gap : 0.;
   checkPressures( couplingScheme, pressure, 1.E-8, "rate" );

   tribol::finalize();
}

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

#ifdef TRIBOL_USE_UMPIRE
  umpire::ResourceManager::getInstance();  // initialize umpire's ResouceManager
#endif

  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}
