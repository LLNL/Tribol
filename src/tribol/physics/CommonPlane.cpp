// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Tribol Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (MIT)

#include "CommonPlane.hpp"

#include "tribol/mesh/MethodCouplingData.hpp"
#include "tribol/mesh/InterfacePairs.hpp"
#include "tribol/mesh/CouplingScheme.hpp"
#include "tribol/geom/ContactPlane.hpp"
#include "tribol/geom/ContactPlaneManager.hpp"
#include "tribol/geom/GeomUtilities.hpp"
#include "tribol/common/Parameters.hpp"
#include "tribol/integ/Integration.hpp"
#include "tribol/integ/FE.hpp"
#include "tribol/utils/ContactPlaneOutput.hpp"
#include "tribol/utils/Math.hpp"

namespace tribol
{

RealT ComputePenaltyStiffnessPerArea( const RealT K1_over_t1,
                                     const RealT K2_over_t2 )
{
   // compute face-pair specific penalty stiffness per unit area.
   // Note: This assumes that each face has a spring stiffness 
   // equal to that side's material Bulk modulus, K, over the 
   // thickness of the volume element to which that face belongs, 
   // times the overlap area. That is, K1_over_t1 * A and K2_over_t2 * A. We 
   // then assume the two springs are in series and compute an 
   // equivalent spring stiffness as, 
   // k_eq = A*(K1_over_t1)*(K2_over_t2) / ((K1_over_t1)+(K2_over_t2). 
   // Note, the host code registers each face's (K/t) as a penalty scale.
   //
   // UNITS: we multiply k_eq above by the overlap area A, to get a 
   // stiffness per unit area. This will make the force calculations 
   // commensurate with the previous calculations using only the 
   // constant registered penalty scale.

   return K1_over_t1 * K2_over_t2 / (K1_over_t1 + K2_over_t2);

} // end ComputePenaltyStiffnessPerArea

//------------------------------------------------------------------------------
RealT ComputeGapRatePressure( ContactPlaneManager& cpMgr, 
                             int cpID, IndexT mesh_id1, IndexT mesh_id2, 
                             int fId1, int fId2, RealT element_penalty,
                             int dim, RatePenaltyCalculation rate_calc )
{
   MeshManager& meshManager = MeshManager::getInstance();
   MeshData& m1 = meshManager.at( mesh_id1 );
   MeshData& m2 = meshManager.at( mesh_id2 );

   // compute the correct rate_penalty
   RealT rate_penalty = 0.;
   switch (rate_calc)
   {
      case NO_RATE_PENALTY:
      {
         return 0.;
      }
      case RATE_CONSTANT:
      {
         rate_penalty = 0.5 * (m1.getElementData().m_rate_penalty_stiffness + 
                               m2.getElementData().m_rate_penalty_stiffness);
         break;
      }
      case RATE_PERCENT:
      {
         rate_penalty = element_penalty * 0.5 * 
                        (m1.getElementData().m_rate_percent_stiffness + 
                         m2.getElementData().m_rate_percent_stiffness);
         break;
      }
      default:
         // no-op, quiet compiler
         break;
   } // end switch on rate_calc

   // compute the velocity gap and pressure contribution
   int numNodesPerCell1 = m1.numberOfNodesPerElement();
   int numNodesPerCell2 = m2.numberOfNodesPerElement();

   RealT x1[dim * numNodesPerCell1];
   RealT v1[dim * numNodesPerCell1];
   m1.getFaceCoords( fId1, &x1[0] );
   m1.getFaceNodalVelocities( fId1, &v1[0] );

   RealT x2[dim * numNodesPerCell2];
   RealT v2[dim * numNodesPerCell2];
   m2.getFaceCoords( fId2, &x2[0] );
   m2.getFaceNodalVelocities( fId2, &v2[0] );

   //////////////////////////////////////////////////////////
   // compute velocity Galerkin approximation at projected // 
   // overlap centroid                                     //
   //////////////////////////////////////////////////////////
   RealT vel_f1[dim];
   RealT vel_f2[dim];
   initRealArray( &vel_f1[0], dim, 0. );
   initRealArray( &vel_f2[0], dim, 0. );

   // interpolate nodal velocity at overlap centroid as projected 
   // onto face 1
   RealT cXf1 = cpMgr.m_cXf1[cpID];
   RealT cYf1 = cpMgr.m_cYf1[cpID];
   RealT cZf1 = (dim == 3) ? cpMgr.m_cZf1[cpID] : 0.;
   GalerkinEval( &x1[0], cXf1, cYf1, cZf1,
                 LINEAR, PHYSICAL, dim, dim, 
                 &v1[0], &vel_f1[0] );

   // interpolate nodal velocity at overlap centroid as projected 
   // onto face 2
   RealT cXf2 = cpMgr.m_cXf2[cpID];
   RealT cYf2 = cpMgr.m_cYf2[cpID];
   RealT cZf2 = (dim == 3) ? cpMgr.m_cZf2[cpID] : 0.;
   GalerkinEval( &x2[0], cXf2, cYf2, cZf2,
                 LINEAR, PHYSICAL, dim, dim, 
                 &v2[0], &vel_f2[0] );

   // compute velocity gap vector
   RealT velGap[dim];
   velGap[0] = vel_f1[0] - vel_f2[0];
   velGap[1] = vel_f1[1] - vel_f2[1];
   if (dim == 3)
   {
      velGap[2] = vel_f1[2] - vel_f2[2];
   }

   // compute velocity gap scalar
   cpMgr.m_velGap[cpID] = 0.;
   cpMgr.m_velGap[cpID] += velGap[0] * cpMgr.m_nX[cpID];
   cpMgr.m_velGap[cpID] += velGap[1] * cpMgr.m_nY[cpID];
   if (dim == 3)
   {
      cpMgr.m_velGap[cpID] += velGap[2] * cpMgr.m_nZ[cpID];
   }

   // check the gap rate sense. 
   // (v1-v2) * \nu < 0 : velocities lead to more interpenetration;
   // note, \nu is in direction of face_2 outward unit normal
   // TODO consider a velocity gap tolerance. Checking this against 
   // 0. actually smoothed out contact behavior in contact problem 1
   // for certain percent rate penalties.
   if (cpMgr.m_velGap[cpID] <= 0.) // TODO do we want = or just <?
   {
      cpMgr.m_ratePressure[ cpID ] = cpMgr.m_velGap[ cpID ] * rate_penalty;
      return cpMgr.m_ratePressure[ cpID ];
   } // end if-check on velocity gap

   return 0.;

} // end ComputeGapRatePressure()

//------------------------------------------------------------------------------
template< >
int ApplyNormal< COMMON_PLANE, PENALTY >( CouplingScheme const * cs )
{
   InterfacePairs const * const pairs = cs->getInterfacePairs();
   const IndexT numPairs = pairs->getNumPairs();

   MeshManager& meshManager = MeshManager::getInstance();
   ContactPlaneManager& cpManager = ContactPlaneManager::getInstance();
   parameters_t& parameters = parameters_t::getInstance();
   int const dim = parameters.dimension;

   LoggingLevel logLevel = cs->getLoggingLevel(); 

   ////////////////////////////////
   // Grab pointers to mesh data //
   ////////////////////////////////
   const IndexT mesh_id1 = cs->getMeshId1();
   const IndexT mesh_id2 = cs->getMeshId2();

   MeshData& mesh1 = meshManager.at( mesh_id1 );
   MeshData& mesh2 = meshManager.at( mesh_id2 );
   const IndexT numNodesPerFace = mesh1.numberOfNodesPerElement();

   auto response1_x = mesh1.getResponse()[0].data();
   auto response1_y = mesh1.getResponse()[1].data();
   RealT* response1_z = nullptr;
   if (mesh1.dimension() == 3)
   {
      response1_z = mesh1.getResponse()[2].data();
   }
   const IndexT * const nodalConnectivity1 = mesh1.getConnectivity().data();

   auto response2_x = mesh2.getResponse()[0].data();
   auto response2_y = mesh2.getResponse()[1].data();
   RealT* response2_z = nullptr;
   if (mesh2.dimension() == 3)
   {
      response2_z = mesh2.getResponse()[2].data();
   }
   const IndexT * nodalConnectivity2 = mesh2.getConnectivity().data();


   ///////////////////////////////
   // loop over interface pairs //
   ///////////////////////////////
   int cpID = 0;
   int err = 0;
   bool neg_thickness {false};
   for (IndexT kp = 0; kp < numPairs; ++kp)
   {
      InterfacePair pair = pairs->getInterfacePair(kp);

      if (!pair.isContactCandidate) 
      {
         continue;
      }

      // get pair indices
      IndexT index1 = pair.pairIndex1;
      IndexT index2 = pair.pairIndex2;

      RealT gap = cpManager.m_gap[ cpID ];
      RealT A = cpManager.m_area[ cpID ]; // face-pair overlap area

      // don't proceed for gaps that don't violate the constraints. This check 
      // allows for numerically zero interpenetration.
      RealT gap_tol = cs->getGapTol( index1, index2 );

      if ( gap > gap_tol )
      {
         // We are here if we have a pair that passes ALL geometric 
         // filter checks, BUT does not actually violate this method's 
         // gap constraint. There will be data in the contact plane 
         // manager so we MUST increment the counter.
         cpManager.m_inContact[ cpID ] = false;
         ++cpID; 
         continue;
      }

      // debug force sums
      RealT dbg_sum_force1 {0.};
      RealT dbg_sum_force2 {0.};

      /////////////////////////////////////////////
      // kinematic penalty stiffness calculation //
      /////////////////////////////////////////////
      RealT penalty_stiff_per_area {0.};
      const EnforcementOptions& enforcement_options = const_cast<EnforcementOptions&>(cs->getEnforcementOptions());
      const PenaltyEnforcementOptions& pen_enfrc_options = enforcement_options.penalty_options;
      RealT pen_scale1 = mesh1.getElementData().m_penalty_scale;
      RealT pen_scale2 = mesh2.getElementData().m_penalty_scale;
      switch (pen_enfrc_options.kinematic_calculation)
      {
         case KINEMATIC_CONSTANT: 
         {
            // pre-multiply each spring stiffness by each mesh's penalty scale
            auto stiffness1 = pen_scale1 * mesh1.getElementData().m_penalty_stiffness;
            auto stiffness2 = pen_scale2 * mesh2.getElementData().m_penalty_stiffness;
            // compute the equivalent contact penalty spring stiffness per area
            penalty_stiff_per_area  = ComputePenaltyStiffnessPerArea( stiffness1, stiffness2 );
            break;
         }
         case KINEMATIC_ELEMENT:
         {
            // add tiny_length to element thickness to avoid division by zero
            auto t1 = mesh1.getElementData().m_thickness[ index1 ] + pen_enfrc_options.tiny_length;
            auto t2 = mesh2.getElementData().m_thickness[ index2 ] + pen_enfrc_options.tiny_length;

            if (t1 < 0. || t2 < 0.)
            {
               neg_thickness = true;
               err = 1;
            }

            // compute each element spring stiffness. Pre-multiply the material modulus 
            // (i.e. material stiffness) by each mesh's penalty scale
            auto stiffness1 = pen_scale1 * mesh1.getElementData().m_mat_mod[ index1 ] / t1;
            auto stiffness2 = pen_scale2 * mesh2.getElementData().m_mat_mod[ index2 ] / t2;
            // compute the equivalent contact penalty spring stiffness per area
            penalty_stiff_per_area = ComputePenaltyStiffnessPerArea( stiffness1, stiffness2 );
            break;
         }
         default:
            // no-op, quiet compiler
            break;
      } // end switch on kinematic penalty calculation option

      ////////////////////////////////////////////////////
      // Compute contact pressure(s) on current overlap // 
      ////////////////////////////////////////////////////

      // compute total pressure based on constraint type
      RealT totalPressure = 0.;
      cpManager.m_pressure[ cpID ] = gap * penalty_stiff_per_area; // kinematic contribution
      switch(pen_enfrc_options.constraint_type)
      {
         case KINEMATIC_AND_RATE:
         {
            // kinematic contribution
            totalPressure += cpManager.m_pressure[ cpID ];
            // add gap-rate contribution
            totalPressure += 
               ComputeGapRatePressure( cpManager, cpID, mesh_id1, mesh_id2, 
                                       index1, index2, penalty_stiff_per_area, dim,
                                       pen_enfrc_options.rate_calculation );
            break;
         }
         case KINEMATIC:
            // kinematic gap pressure contribution  only
            totalPressure += cpManager.m_pressure[ cpID ];
            break;
         default:
            // no-op
            break;
      } // end switch on registered penalty enforcement option

      // debug prints. Comment out for now, but keep for future common plane 
      // debugging
//         SLIC_DEBUG("gap: " << gap);
//         SLIC_DEBUG("area: " << A);
//         SLIC_DEBUG("penalty stiffness: " << penalty_stiff_per_area);
//         SLIC_DEBUG("pressure: " << cpManager.m_pressure[ cpID ]);

      ///////////////////////////////////////////
      // create surface contact element struct //
      ///////////////////////////////////////////

      // construct array of nodal coordinates
      RealT xf1[ dim * numNodesPerFace ];
      RealT xf2[ dim * numNodesPerFace ];
      RealT xVert[dim * cpManager.m_numPolyVert[cpID]];  

      initRealArray( &xf1[0], dim*numNodesPerFace, 0. );
      initRealArray( &xf2[0], dim*numNodesPerFace, 0. );
      initRealArray( &xVert[0], dim*cpManager.m_numPolyVert[cpID], 0. );

//      // get projected face coordinates
//      cpManager.getProjectedFaceCoords( cpID, 0, &xf1[0] ); // face 0 = first face
//      cpManager.getProjectedFaceCoords( cpID, 1, &xf2[0] ); // face 1 = second face

      // get current configuration, physical coordinates of each face
      mesh1.getFaceCoords( index1, &xf1[0] );
      mesh2.getFaceCoords( index2, &xf2[0] );

      // construct array of polygon overlap vertex coordinates
      cpManager.getContactPlaneOverlapVerts( cpID, cpManager.m_numPolyVert[cpID], 
                                             &xVert[0] );

      // instantiate surface contact element struct. Note, this is done with current 
      // configuration face coordinates (i.e. NOT on the contact plane) and overlap 
      // coordinates ON the contact plane. The surface contact element does not need 
      // to be used this way, but the developer should do the book-keeping.
      SurfaceContactElem cntctElem( dim, &xf1[0], &xf2[0], &xVert[0],
                                    numNodesPerFace, 
                                    cpManager.m_numPolyVert[cpID],
                                    mesh_id1, mesh_id2, index1, index2 );

      // set SurfaceContactElem face normals and overlap normal
      RealT faceNormal1[dim];
      RealT faceNormal2[dim];
      RealT overlapNormal[dim];

      mesh1.getFaceNormal( index1, dim, &faceNormal1[0] );
      mesh2.getFaceNormal( index2, dim, &faceNormal2[0] );
      cpManager.getContactPlaneNormal( cpID, dim, &overlapNormal[0] );

      cntctElem.faceNormal1 = &faceNormal1[0];
      cntctElem.faceNormal2 = &faceNormal2[0];
      cntctElem.overlapNormal = &overlapNormal[0];

      // create arrays to hold nodal residual weak form integral evaluations
      RealT phi1[numNodesPerFace];
      RealT phi2[numNodesPerFace];
      initRealArray( &phi1[0], numNodesPerFace, 0. );
      initRealArray( &phi2[0], numNodesPerFace, 0. );

      ////////////////////////////////////////////////////////////////////////
      // Integration of contact integrals: integral of shape functions over //
      // contact overlap patch                                              //
      ////////////////////////////////////////////////////////////////////////
      EvalWeakFormIntegral< COMMON_PLANE, SINGLE_POINT >
                          ( cntctElem, &phi1[0], &phi2[0] );

      ///////////////////////////////////////////////////////////////////////
      // Computation of full contact nodal force contributions             //
      // (i.e. premultiplication of contact integrals by normal component, //
      //  contact pressure, and overlap area)                              //
      ///////////////////////////////////////////////////////////////////////

      RealT phi_sum_1 = 0.;
      RealT phi_sum_2 = 0.;

      // compute contact force (spring force)
      RealT contact_force = totalPressure * A;
      RealT force_x = overlapNormal[0] * contact_force;
      RealT force_y = overlapNormal[1] * contact_force;
      RealT force_z = 0.;
      if (dim == 3)
      {
         force_z = overlapNormal[2] * contact_force;
      }

      //////////////////////////////////////////////////////
      // loop over nodes and compute contact nodal forces //
      //////////////////////////////////////////////////////
      for( IndexT a=0 ; a < numNodesPerFace ; ++a )
      {

        IndexT node0 = nodalConnectivity1[ index1*numNodesPerFace + a ];
        IndexT node1 = nodalConnectivity2[ index2*numNodesPerFace + a ];

        if (logLevel == TRIBOL_DEBUG)
        {
           phi_sum_1 += phi1[a];
           phi_sum_2 += phi2[a];
        }
 
        const RealT nodal_force_x1 = force_x * phi1[a];
        const RealT nodal_force_y1 = force_y * phi1[a];
        const RealT nodal_force_z1 = force_z * phi1[a];

        const RealT nodal_force_x2 = force_x * phi2[a];
        const RealT nodal_force_y2 = force_y * phi2[a];
        const RealT nodal_force_z2 = force_z * phi2[a];

        if (logLevel == TRIBOL_DEBUG)
        {
           dbg_sum_force1 += magnitude( nodal_force_x1, 
                                        nodal_force_y1, 
                                        nodal_force_z1 );
           dbg_sum_force2 += magnitude( nodal_force_x2,
                                        nodal_force_y2,
                                        nodal_force_z2 );
        }

        // accumulate contributions in host code's registered nodal force arrays
        response1_x[ node0 ] -= nodal_force_x1;
        response2_x[ node1 ] += nodal_force_x2;

        response1_y[ node0 ] -= nodal_force_y1;
        response2_y[ node1 ] += nodal_force_y2;

        // there is no z component for 2D
        if (mesh1.dimension() == 3)
        {
           response1_z[ node0 ] -= nodal_force_z1;
           response2_z[ node1 ] += nodal_force_z2;
        }
      } // end for loop over face nodes

      // comment out debug logs; too much output during tests. Keep for easy 
      // debugging if needed
      //SLIC_DEBUG("force sum, side 1, pair " << kp << ": " << -dbg_sum_force1 );
      //SLIC_DEBUG("force sum, side 2, pair " << kp << ": " << dbg_sum_force2 );
      //SLIC_DEBUG("phi 1 sum: " << phi_sum_1 );
      //SLIC_DEBUG("phi 2 sum: " << phi_sum_2 );
    
      // increment contact plane id 
      ++cpID;

   } // end loop over interface pairs
  
   SLIC_DEBUG_IF(neg_thickness, "ApplyNormal<COMMON_PLANE, PENALTY>: negative element thicknesses encountered.");

   return err;

} // end ApplyNormal<COMMON_PLANE, PENALTY>()

}
