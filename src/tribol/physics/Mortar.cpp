// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Tribol Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (MIT)

#include "Mortar.hpp"

#include "tribol/types.hpp"
#include "tribol/mesh/MeshManager.hpp"
#include "tribol/mesh/MethodCouplingData.hpp"
#include "tribol/mesh/InterfacePairs.hpp"
#include "tribol/mesh/CouplingScheme.hpp"
#include "tribol/mesh/CouplingSchemeManager.hpp"
#include "tribol/geom/ContactPlane.hpp"
#include "tribol/geom/ContactPlaneManager.hpp"
#include "tribol/geom/GeomUtilities.hpp"
#include "tribol/common/Parameters.hpp"
#include "tribol/integ/Integration.hpp"
#include "tribol/integ/FE.hpp"
#include "tribol/utils/ContactPlaneOutput.hpp"
#include "tribol/utils/Math.hpp"
#include "tribol/common/logger.hpp"

// Axom includes
#include "axom/slic.hpp"

#include <iostream>

namespace tribol
{

void ComputeMortarWeights( SurfaceContactElem & elem )
{
   // instantiate integration object
   IntegPts integ;

   // Debug: leave code in for now to call Gauss quadrature on triangle rule
   GaussPolyIntTri( elem, integ, 3 );

   // call Taylor-Wingate-Bos integation rule. NOTE: this is not 
   // working. The correct gaps are not being computed.
//   TWBPolyInt( elem, integ, 3 );

   // get individual arrays of coordinates for each face
   real x1[elem.numFaceVert];
   real y1[elem.numFaceVert];
   real z1[elem.numFaceVert];
   real x2[elem.numFaceVert];
   real y2[elem.numFaceVert];
   real z2[elem.numFaceVert];

   for (int i=0; i<elem.numFaceVert; ++i)
   {
      x1[i] = elem.faceCoords1[elem.dim*i];
      y1[i] = elem.faceCoords1[elem.dim*i+1];
      z1[i] = elem.faceCoords1[elem.dim*i+2];
      x2[i] = elem.faceCoords2[elem.dim*i];
      y2[i] = elem.faceCoords2[elem.dim*i+1];
      z2[i] = elem.faceCoords2[elem.dim*i+2];
   }

   // allocate mortar weights array on SurfaceContactElem object. This routine 
   // also initializes the array
   elem.allocateMortarWts();

   real phiSlaveA, phiSlaveB, phiMasterA;

   // loop over number of nodes on the slave or master depending on whether forming 
   // slave/slave or master/slave weights
  
   for (int a=0; a<elem.numFaceVert; ++a)
   {
      // loop over number of nodes on slave side
      for (int b=0; b<elem.numFaceVert; ++b)
      {
         // set slave/slave and master/slave ids...Don't change these ids
         int slaveSlaveId = elem.numFaceVert * a + b;
         int masterSlaveId = elem.numFaceVert * elem.numFaceVert + elem.numFaceVert * a + b;

         // loop over number of integration points
         for (int ip=0; ip<integ.numIPs; ++ip)
         {
            // The integration method for computing weights uses 
            // the inverse isoparametric mapping of a current configuration 
            // integration point (as projected onto the current configuration 
            // face) to obtain a (xi,eta) coordinate pair in parent space 
            // for the evaluation of Lagrange shape functions
            real xp[3] = { integ.xy[elem.dim*ip], integ.xy[elem.dim*ip+1], integ.xy[elem.dim*ip+2] };
            real xi[2] = { 0., 0. };

            InvIso( xp, x1, y1, z1, elem.numFaceVert, xi );
            LinIsoQuadShapeFunc( xi[0], xi[1], a, phiMasterA );

            InvIso( xp, x2, y2, z2, elem.numFaceVert, xi );
            LinIsoQuadShapeFunc( xi[0], xi[1], a, phiSlaveA );
            LinIsoQuadShapeFunc( xi[0], xi[1], b, phiSlaveB );

            #ifdef TRIBOL_DEBUG_LOG
               if (slaveSlaveId > elem.numWts || masterSlaveId > elem.numWts)
               {
                  TRIBOL_ERROR("ComputeMortarWts: integer ids for weights exceed elem.numWts");
               }
            #endif /* TRIBOL_DEBUG_LOG */

            // compute slave/slave mortar weight
            elem.mortarWts[ slaveSlaveId ]  += integ.wts[ip] * phiSlaveA * phiSlaveB;


            // compute master/slave mortar weight
            elem.mortarWts[ masterSlaveId ] += integ.wts[ip] * phiMasterA * phiSlaveB;

         } // end loop over integration points

      } // end loop over nodes on side 2

   } // end loop over nodes on side 1
   
} // end ComputeMortarWeights()

//------------------------------------------------------------------------------
template< >
void ComputeNodalGap< SINGLE_MORTAR >( SurfaceContactElem & elem )
{
   // check to make sure mortar weights have been computed locally 
   // for the SurfaceContactElem object
   #ifdef TRIBOL_DEBUG_LOG
      if (elem.mortarWts == nullptr)
      {
         TRIBOL_ERROR("ComputeNodalGap< SINGLE_MORTAR >: compute local weights on input struct first.");
      }
   #endif /* TRIBOL_DEBUG_LOG */

   // get mesh instance to store gaps on mesh data object
   MeshManager& meshManager = MeshManager::getInstance();
   MeshData& slaveMesh = meshManager.GetMeshInstance( elem.meshId2 );
   IndexType const * const slaveConn = slaveMesh.m_connectivity;

   // will populate local gaps on slave face on slave mesh data object
   if (slaveMesh.m_nodalFields.m_node_gap == nullptr)
   {
      SLIC_ERROR("ComputeNodalGap< SINGLE_MORTAR >: allocate gaps on mesh data object."); 
   }

   if (slaveMesh.m_node_nX == nullptr || slaveMesh.m_node_nY == nullptr)
   {
      SLIC_ERROR("ComputeNodalGap< SINGLE_MORTAR >: allocate and compute nodal normals on mesh data object.");   
   }

   // compute gap contributions associated with face 2 on the SurfaceContactElem 
   // (i.e. slave surface)

   // loop over nodes on slave side
   for (int a=0; a<elem.numFaceVert; ++a)
   {
      // initialize gap1 and gap2 terms
      real g1 = 0.;
      real g2 = 0.;

      // get global slave node number from connectivity
      real nrml_a[elem.dim];
      int glbId = slaveConn[ elem.numFaceVert * elem.faceId2 + a ];
      nrml_a[0] = slaveMesh.m_node_nX[ glbId ];
      nrml_a[1] = slaveMesh.m_node_nY[ glbId ];
      if (elem.dim == 3 )
      {
         nrml_a[2] = slaveMesh.m_node_nZ[ glbId ];
      }

      // sum contributions from both sides
      for (int b=0; b<elem.numFaceVert; ++b)
      {
         // compute slave-master and slave-slave ids. Note, n_ab is 
         // the stored mortar weight. For master-slave mortar weights, 
         // a = master node and b = slave node, BUT FOR THE GAP COMPUTATION,
         // THE SUM OF MORTAR WEIGHTS IS ACTUALLY OVER SHAPE FUNCTIONS 
         // DEFINED AT NODE "b", SO WE NEED TO USE (n_ab)^T.
         real nab_1 = elem.getSlaveMasterWt( a, b ); // slave-master weight
         real nab_2 = elem.getSlaveSlaveWt( a, b ); // slave-slave weight

         g1 += dotProd( &nrml_a[0], &elem.faceCoords1[ elem.dim * b ], elem.dim ) *
               nab_1;
         g2 += dotProd( &nrml_a[0], &elem.faceCoords2[ elem.dim * b ], elem.dim ) * 
               nab_2;
      }

      // store local gap
      slaveMesh.m_nodalFields.m_node_gap[ glbId ] += (g1-g2);

   } // end a-loop over slave nodes

} // end ComputeNodalGap<>()

//------------------------------------------------------------------------------
void ComputeSingleMortarGaps( CouplingScheme const * cs )
{
   InterfacePairs const * const pairs = cs->getInterfacePairs();
   IndexType const numPairs = pairs->getNumPairs();

   MeshManager& meshManager = MeshManager::getInstance();
   ContactPlaneManager& cpManager = ContactPlaneManager::getInstance();
   parameters_t& parameters = parameters_t::getInstance();
   integer const dim = parameters.dimension;
   IndexType const numNodesPerFace = (dim == 3) ? 4 : 2;

   ////////////////////////////////////////////////////////////////////////
   //
   // Grab pointers to mesh data
   //
   ////////////////////////////////////////////////////////////////////////
   IndexType const masterId = cs->getMeshId1();
   IndexType const slaveId = cs->getMeshId2();

   MeshData& masterMesh = meshManager.GetMeshInstance( masterId );
   MeshData& slaveMesh = meshManager.GetMeshInstance( slaveId );

   real const * const x1 = masterMesh.m_positionX;
   real const * const y1 = masterMesh.m_positionY; 
   real const * const z1 = masterMesh.m_positionZ; 
   IndexType const * const masterConn= masterMesh.m_connectivity;

   real const * const x2 = slaveMesh.m_positionX; 
   real const * const y2 = slaveMesh.m_positionY;
   real const * const z2 = slaveMesh.m_positionZ;
   IndexType const * slaveConn = slaveMesh.m_connectivity;

   // compute nodal normals (do this outside the element loop)
   slaveMesh.computeNodalNormals( dim );

   // declare local variables to hold face nodal coordinates
   // and overlap vertex coordinates
   IndexType size = dim * numNodesPerFace;
   real masterX[ size ];
   real slaveX[ size ];

   // projected coords
   real masterX_bar[ size ];
   real slaveX_bar[ size ];
   real* overlapX;

   ////////////////////////////////////////////////////////////////
   // compute slave gaps to determine active set of contact dofs //
   ////////////////////////////////////////////////////////////////
   int cpID = 0;
   for (IndexType kp = 0; kp < numPairs; ++kp)
   {
      InterfacePair pair = pairs->getInterfacePair(kp);

      if (!pair.inContact)
      {
         continue;
      }

      // get pair indices
      IndexType index1 = pair.pairIndex1;
      IndexType index2 = pair.pairIndex2;

      // populate the current configuration nodal coordinates for the 
      // two faces
      for (int i=0; i<numNodesPerFace; ++i)
      {
         int id = dim * i;
         IndexType master_id = masterConn[ numNodesPerFace * index1 + i ];
         IndexType slave_id  = slaveConn[ numNodesPerFace * index2 + i ];

         masterX[ id ]   = x1[ master_id ];
         masterX[ id+1 ] = y1[ master_id ];
         masterX[ id+2 ] = z1[ master_id ];
         slaveX[ id ]   = x2[ slave_id ];
         slaveX[ id+1 ] = y2[ slave_id ];
         slaveX[ id+2 ] = z2[ slave_id ];
    
         // arrays for projected coords
         masterX_bar[ id ]   = x1[ master_id ];
         masterX_bar[ id+1 ] = y1[ master_id ];
         masterX_bar[ id+2 ] = z1[ master_id ];
         slaveX_bar[ id ]   = x2[ slave_id ];
         slaveX_bar[ id+1 ] = y2[ slave_id ];
         slaveX_bar[ id+2 ] = z2[ slave_id ];
      }

      overlapX = new real[ dim * cpManager.m_numPolyVert[ cpID ]];
      initRealArray( overlapX, dim*cpManager.m_numPolyVert[cpID], 0. );

      // get projected face coordinates
      cpManager.getProjectedFaceCoords( cpID, 0, &masterX_bar[0] ); // face 0 = first face
      cpManager.getProjectedFaceCoords( cpID, 1, &slaveX_bar[0] ); // face 1 = second face

      // construct array of polygon overlap vertex coordinates
      cpManager.getContactPlaneOverlapVerts( cpID, cpManager.m_numPolyVert[cpID], &overlapX[0] );

      // instantiate contact surface element for purposes of computing 
      // mortar weights. Note, this uses projected face coords
      SurfaceContactElem elem( dim, &masterX_bar[0], &slaveX_bar[0], 
                               &overlapX[0],
                               numNodesPerFace, 
                               cpManager.m_numPolyVert[cpID],
                               masterId, slaveId, index1, index2 );

      // compute the mortar weights to be stored on the surface 
      // contact element struct. This must be done prior to computing nodal gaps
      elem.overlapArea = cpManager.m_area[ cpID ];
      ComputeMortarWeights( elem );

      // compute mortar gaps. Note, we have to now use current configuration
      // nodal coordinates on the contact element
      elem.faceCoords1 = &masterX[0];
      elem.faceCoords2 = &slaveX[0];

      ComputeNodalGap< SINGLE_MORTAR >( elem );

      // TODO: fix this to register the actual number of active slave gaps.
      // This is not the appropriate data structure to put this information in 
      // as the SurfaceContactElem goes out of scope when we exit the loop.
      // HAVE TO set the number of active constraints. For now set to 
      // all slave face nodes.
      elem.numActiveGaps = numNodesPerFace;

      ++cpID;

      delete [] overlapX;

   } // end loop over pairs to compute nodal gaps

   slaveMesh.m_nodalFields.m_isGapComputed = true;

} // end ComputeSingleMortarGaps()

//------------------------------------------------------------------------------
template< >
int ApplyNormal< SINGLE_MORTAR, LAGRANGE_MULTIPLIER >( CouplingScheme const * cs )
{
   InterfacePairs const * const pairs = cs->getInterfacePairs();
   IndexType const numPairs = pairs->getNumPairs();

   MeshManager& meshManager = MeshManager::getInstance();
   ContactPlaneManager& cpManager = ContactPlaneManager::getInstance();
   parameters_t& parameters = parameters_t::getInstance();
   integer const dim = parameters.dimension;
   IndexType const numNodesPerFace = (dim == 3) ? 4 : 2;

   ////////////////////////////////
   //                            //
   // Grab pointers to mesh data //
   //                            //
   ////////////////////////////////
   IndexType const masterId = cs->getMeshId1();
   IndexType const slaveId = cs->getMeshId2();

   MeshData& masterMesh = meshManager.GetMeshInstance( masterId );
   MeshData& slaveMesh = meshManager.GetMeshInstance( slaveId );

   real const * const x1 = masterMesh.m_positionX;
   real const * const y1 = masterMesh.m_positionY; 
   real const * const z1 = masterMesh.m_positionZ; 
   real * const fx1 = masterMesh.m_forceX;
   real * const fy1 = masterMesh.m_forceY; 
   real * const fz1 = masterMesh.m_forceZ; 
   IndexType const * const masterConn= masterMesh.m_connectivity;

   real const * const x2 = slaveMesh.m_positionX; 
   real const * const y2 = slaveMesh.m_positionY;
   real const * const z2 = slaveMesh.m_positionZ;
   real * const fx2 = slaveMesh.m_forceX; 
   real * const fy2 = slaveMesh.m_forceY;
   real * const fz2 = slaveMesh.m_forceZ;
   IndexType const * slaveConn = slaveMesh.m_connectivity;

   // projected coords
   real masterX_bar[ dim * numNodesPerFace ];
   real slaveX_bar[ dim * numNodesPerFace ];
   real* overlapX; // [dim * cpManager.m_numPolyVert[cpID]];  

   ////////////////////////////////
   //                            //
   // compute single mortar gaps //
   //                            //
   ////////////////////////////////
   ComputeSingleMortarGaps( cs );

   int numTotalNodes = cs->getNumTotalNodes();
   int numRows = dim * numTotalNodes + numTotalNodes;
   const EnforcementOptions& enforcement_options = const_cast<EnforcementOptions&>(cs->getEnforcementOptions());
   const LagrangeMultiplierImplicitOptions& lm_options  = enforcement_options.lm_implicit_options;
   if ( lm_options.sparse_mode == SparseMode::MFEM_ELEMENT_DENSE )
   {
      static_cast<MortarData*>( cs->getMethodData() )->reserveBlockJ( 
         {BlockSpace::MASTER, BlockSpace::SLAVE, BlockSpace::LAGRANGE_MULTIPLIER},
         numPairs
      );
   }
   else if ( lm_options.sparse_mode == SparseMode::MFEM_INDEX_SET || 
             lm_options.sparse_mode == SparseMode::MFEM_LINKED_LIST )
   {
      static_cast<MortarData*>( cs->getMethodData() )->allocateMfemSparseMatrix( numRows );
   }
   else
   {
      SLIC_ERROR("Unsupported Jacobian storage method.");
   }

   ////////////////////////////////////////////////////////////////
   //                                                            //
   // compute equilibrium residual and/or Jacobian contributions //
   //                                                            //
   ////////////////////////////////////////////////////////////////
   int cpID = 0;
   for (IndexType kp = 0; kp < numPairs; ++kp)
   {
      InterfacePair pair = pairs->getInterfacePair(kp);

      if (!pair.inContact)
      {
         continue;
      }

      // get pair indices
      IndexType index1 = pair.pairIndex1;
      IndexType index2 = pair.pairIndex2;

      // populate the current configuration nodal coordinates for the 
      // two faces
      for (int i=0; i<numNodesPerFace; ++i)
      {
         int id = dim * i;
         IndexType master_id = masterConn[ numNodesPerFace * index1 + i];
         IndexType slave_id = slaveConn[ numNodesPerFace * index2 + i ]; 
         // arrays for projected coords
         masterX_bar[ id ]   = x1[ master_id ];
         masterX_bar[ id+1 ] = y1[ master_id ];
         masterX_bar[ id+2 ] = z1[ master_id ];
         slaveX_bar[ id ]   = x2[ slave_id ];
         slaveX_bar[ id+1 ] = y2[ slave_id ];
         slaveX_bar[ id+2 ] = z2[ slave_id ];
      }

      overlapX = new real[ dim * cpManager.m_numPolyVert[ cpID ]];
      initRealArray( overlapX, dim*cpManager.m_numPolyVert[cpID], 0. );

      // get projected face coordinates
      cpManager.getProjectedFaceCoords( cpID, 0, &masterX_bar[0] ); // face 0 = first face
      cpManager.getProjectedFaceCoords( cpID, 1, &slaveX_bar[0] ); // face 1 = second face

      // construct array of polygon overlap vertex coordinates
      cpManager.getContactPlaneOverlapVerts( cpID, cpManager.m_numPolyVert[cpID], &overlapX[0] );

      // instantiate contact surface element for purposes of computing 
      // mortar weights. Note, this uses projected face coords
      SurfaceContactElem elem( dim, &masterX_bar[0], &slaveX_bar[0], 
                               &overlapX[0],
                               numNodesPerFace, 
                               cpManager.m_numPolyVert[cpID],
                               masterId, slaveId, index1, index2 );

      //////////////////////////////////
      // compute equilibrium residual //
      //////////////////////////////////

      // compute mortar weight
      elem.overlapArea = cpManager.m_area[ cpID ];
      ComputeMortarWeights( elem );

      // TODO fix this. This may not be required.
      // HAVE TO set the number of active constraints. For now set to 
      // all slave face nodes.
      elem.numActiveGaps = numNodesPerFace;

      // loop over face nodes (BOTH MASTER and SLAVE 
      // contributions)
      for (int a=0; a<numNodesPerFace; ++a)
      {
         int masterIdA = masterConn[ index1 * numNodesPerFace + a];
         int slaveIdA = slaveConn[ index2 * numNodesPerFace + a ];

         // inner loop over SLAVE nodes
         for (int b=0; b<numNodesPerFace; ++b)
         {
            int slaveIdB = slaveConn[ index2 * numNodesPerFace + b ];

            // We include all slave nodes even if nodal gap is in separation. 
            // NOTE: Per testing, we include ALL slave nodes 
            // in the computation after the geometric filtering and judge contact 
            // activity based on the gap AND the pressure solution

            real forceX = slaveMesh.m_nodalFields.m_node_pressure[ slaveIdB ] * 
                          slaveMesh.m_node_nX[ slaveIdB ];
            real forceY = slaveMesh.m_nodalFields.m_node_pressure[ slaveIdB ] * 
                          slaveMesh.m_node_nY[ slaveIdB ];
            real forceZ = slaveMesh.m_nodalFields.m_node_pressure[ slaveIdB ] * 
                          slaveMesh.m_node_nZ[ slaveIdB ];

            // contact nodal force is the interpolated force using mortar 
            // weights n_ab, where "a" is master or slave node and "b" is 
            // slave node.
            fx1[ masterIdA ] += forceX * elem.getMasterSlaveWt( a, b );
            fy1[ masterIdA ] += forceY * elem.getMasterSlaveWt( a, b ); 
            fz1[ masterIdA ] += forceZ * elem.getMasterSlaveWt( a, b ); 

            fx2[ slaveIdA ]  -= forceX * elem.getSlaveSlaveWt( a, b );
            fy2[ slaveIdA ]  -= forceY * elem.getSlaveSlaveWt( a, b );
            fz2[ slaveIdA ]  -= forceZ * elem.getSlaveSlaveWt( a, b );

         } // end inner loop over slave nodes

      } // end outer loop over slave and master nodes

      //////////////////////////////////////////////////////////
      // compute tangent stiffness contributions if requested //
      //////////////////////////////////////////////////////////
      if ( lm_options.eval_mode == ImplicitEvalMode::MORTAR_RESIDUAL_JACOBIAN ||
           lm_options.eval_mode == ImplicitEvalMode::MORTAR_JACOBIAN )
      {
         ComputeSingleMortarJacobian( elem );
         if ( lm_options.sparse_mode == SparseMode::MFEM_ELEMENT_DENSE )
         {
            static_cast<MortarData*>( cs->getMethodData() )->storeElemBlockJ(
               {elem.faceId1, elem.faceId2, elem.faceId2},
               elem.blockJ
            );
         }
         else if ( lm_options.sparse_mode == SparseMode::MFEM_INDEX_SET || 
                   lm_options.sparse_mode == SparseMode::MFEM_LINKED_LIST )
         {
            static_cast<MortarData*>( cs->getMethodData() )->assembleJacobian( elem, lm_options.sparse_mode );
         }
         else
         {
            SLIC_ERROR("Unsupported Jacobian storage method.");
         }
         
      }

      ++cpID;

      delete [] overlapX;

   } // end of loop over interface pairs computing residual/Jacobian contributions

   return 0;

} // end ApplyNormal<>()

//------------------------------------------------------------------------------
template< >
void ComputeResidualJacobian< SINGLE_MORTAR, PRIMAL >( SurfaceContactElem & TRIBOL_UNUSED_PARAM(elem) )
{
   // There is no Jacobian contribution for this block. Be safe and zero out...
   return;
}

//------------------------------------------------------------------------------
template< >
void ComputeResidualJacobian< SINGLE_MORTAR, DUAL >( SurfaceContactElem & elem )
{
   MeshManager& meshManager = MeshManager::getInstance();
   MeshData& slaveMesh = meshManager.GetMeshInstance( elem.meshId2 );
   IndexType const * const slaveConn = slaveMesh.m_connectivity;

   // loop over "a" nodes accumulating sums of master/slave 
   // and slave/slave weights
   for (int a = 0; a<elem.numFaceVert; ++a)
   {
      // single loop over "b" nodes accumulating sums of 
      // master(a)/slave(b) and slave(a)/slave(b) weights
      for (int b = 0; b<elem.numFaceVert; ++b)
      {
         // get global slave node id to index into nodal normals on 
         // slave mesh
         real nrml_b[elem.dim];
         int glbId = slaveConn[ elem.numFaceVert * elem.faceId2 + b ];

         // We assemble ALL slave node contributions, even if gap is in separation.
         // NOTE: Per testing, we compute ALL slave nodes 
         // for faces that have positive areas of overlap after the geometric 
         // filtering and use the gap AND the pressure solution to determine 
         // contact activity

         nrml_b[0] = slaveMesh.m_node_nX[ glbId ];
         nrml_b[1] = slaveMesh.m_node_nY[ glbId ];
         if (elem.dim == 3 )
         {
            nrml_b[2] = slaveMesh.m_node_nZ[ glbId ];
         }

         // get master-slave and slave-slave mortar weights
         real n_master_b = elem.getMasterSlaveWt( a, b ); // master-slave weight
         real n_slave_b  = elem.getSlaveSlaveWt( a, b ); // slave-slave weight, note negative in formulation
         
         // fill Jrp element-pair Jacobian blocks
         // Fill block (0, 2)
         int elem_xdof = elem.getJacobianIndex(SurfaceContactElem::JrpBlock, a, b );
         int dim_offset = elem.getJacobianDimOffset(SurfaceContactElem::JrpBlock);
         elem.blockJ(
            static_cast<axom::IndexType>(BlockSpace::MASTER),
            static_cast<axom::IndexType>(BlockSpace::LAGRANGE_MULTIPLIER)
         ).Data()[ elem_xdof ]                += nrml_b[0] * n_master_b;
         elem.blockJ(
            static_cast<axom::IndexType>(BlockSpace::MASTER),
            static_cast<axom::IndexType>(BlockSpace::LAGRANGE_MULTIPLIER)
         ).Data()[ elem_xdof + dim_offset ]   += nrml_b[1] * n_master_b;
         elem.blockJ(
            static_cast<axom::IndexType>(BlockSpace::MASTER),
            static_cast<axom::IndexType>(BlockSpace::LAGRANGE_MULTIPLIER)
         ).Data()[ elem_xdof + 2*dim_offset ] += nrml_b[2] * n_master_b;

         // Fill block (1, 2)
         elem.blockJ(
            static_cast<axom::IndexType>(BlockSpace::SLAVE),
            static_cast<axom::IndexType>(BlockSpace::LAGRANGE_MULTIPLIER)
         ).Data()[ elem_xdof ]                -= nrml_b[0] * n_slave_b;
         elem.blockJ(
            static_cast<axom::IndexType>(BlockSpace::SLAVE),
            static_cast<axom::IndexType>(BlockSpace::LAGRANGE_MULTIPLIER)
         ).Data()[ elem_xdof + dim_offset ]   -= nrml_b[1] * n_slave_b;
         elem.blockJ(
            static_cast<axom::IndexType>(BlockSpace::SLAVE),
            static_cast<axom::IndexType>(BlockSpace::LAGRANGE_MULTIPLIER)
         ).Data()[ elem_xdof + 2*dim_offset ] -= nrml_b[2] * n_slave_b;

      } // end loop over b nodes

   } // end loop over a nodes

   return;
} // end ComputeResidualJacobian<>()

//------------------------------------------------------------------------------
template< >
void ComputeConstraintJacobian< SINGLE_MORTAR, PRIMAL >( SurfaceContactElem & elem )
{
   MeshManager& meshManager = MeshManager::getInstance();
   MeshData& slaveMesh = meshManager.GetMeshInstance( elem.meshId2 );
   IndexType const * const slaveConn = slaveMesh.m_connectivity;

   // loop over slave nodes for which we are accumulating Jacobian 
   // contributions
   for (int a = 0; a<elem.numFaceVert; ++a)
   {
      // get global slave node id to index into nodal normals on 
      // slave mesh
      real nrml_a[elem.dim];
      int glbId = slaveConn[ elem.numFaceVert * elem.faceId2 + a ];

      // We assemble ALL slave node contributions even if gap is in separation.
      // NOTE: Per Mike Puso we compute ALL slave node 
      // contributions for faces that have positive areas of overlap per the 
      // geometric filtering. Contact activity is judged based on gaps AND 
      // the pressure solution.

      nrml_a[0] = slaveMesh.m_node_nX[ glbId ];
      nrml_a[1] = slaveMesh.m_node_nY[ glbId ];
      if (elem.dim == 3 )
      {
         nrml_a[2] = slaveMesh.m_node_nZ[ glbId ];
      }

      // single loop over "b" nodes accumulating sums of 
      // slave(a)/master(b) and slave(a)/slave(b) weights
      for (int b = 0; b<elem.numFaceVert; ++b)
      {
         // get slave-master and slave-slave mortar weights
         real n_master_a = elem.getSlaveMasterWt( a, b ); // slave-master weight
         real n_slave_a  = elem.getSlaveSlaveWt( a, b ); // slave-slave weight, note negative in formulation

         // fill Jgu element-pair Jacobian blocks
         // Fill block (2, 0)
         int dim_offset = elem.getJacobianDimOffset(SurfaceContactElem::JguBlock);
         int elem_xdof = elem.getJacobianIndex(SurfaceContactElem::JguBlock, a, b );
         elem.blockJ(
            static_cast<axom::IndexType>(BlockSpace::LAGRANGE_MULTIPLIER),
            static_cast<axom::IndexType>(BlockSpace::MASTER)
         ).Data()[ elem_xdof ]                += nrml_a[0] * n_master_a;
         elem.blockJ(
            static_cast<axom::IndexType>(BlockSpace::LAGRANGE_MULTIPLIER),
            static_cast<axom::IndexType>(BlockSpace::MASTER)
         ).Data()[ elem_xdof + dim_offset ]   += nrml_a[1] * n_master_a;
         elem.blockJ(
            static_cast<axom::IndexType>(BlockSpace::LAGRANGE_MULTIPLIER),
            static_cast<axom::IndexType>(BlockSpace::MASTER)
         ).Data()[ elem_xdof + 2*dim_offset ] += nrml_a[2] * n_master_a;

         // Fill block (2, 1)
         elem.blockJ(
            static_cast<axom::IndexType>(BlockSpace::LAGRANGE_MULTIPLIER),
            static_cast<axom::IndexType>(BlockSpace::SLAVE)
         ).Data()[ elem_xdof ]                -= nrml_a[0] * n_slave_a;
         elem.blockJ(
            static_cast<axom::IndexType>(BlockSpace::LAGRANGE_MULTIPLIER),
            static_cast<axom::IndexType>(BlockSpace::SLAVE)
         ).Data()[ elem_xdof + dim_offset ]   -= nrml_a[1] * n_slave_a;
         elem.blockJ(
            static_cast<axom::IndexType>(BlockSpace::LAGRANGE_MULTIPLIER),
            static_cast<axom::IndexType>(BlockSpace::SLAVE)
         ).Data()[ elem_xdof + 2*dim_offset ] -= nrml_a[2] * n_slave_a;

      } // end loop over b nodes

   } // end loop over a nodes

   return;
} // end ComputeConstraintJacobian

//------------------------------------------------------------------------------
template< >
void ComputeConstraintJacobian< SINGLE_MORTAR, DUAL >( SurfaceContactElem& TRIBOL_UNUSED_PARAM(elem) )
{
   // unless we end up solving the complementarity equation, there is 
   // no Jacobian contribtion for this block. Zero out to be safe...
   return;
}

//------------------------------------------------------------------------------
void ComputeSingleMortarJacobian( SurfaceContactElem & elem )
{
   elem.allocateBlockJ( LAGRANGE_MULTIPLIER );

   ComputeResidualJacobian  < SINGLE_MORTAR, PRIMAL >( elem );

   ComputeResidualJacobian  < SINGLE_MORTAR, DUAL   >( elem );

   ComputeConstraintJacobian< SINGLE_MORTAR, PRIMAL >( elem );

   ComputeConstraintJacobian< SINGLE_MORTAR, DUAL   >( elem );

   // Optionally print contact element matrix. Keep commented out here.
   //elem.printBlockJMatrix();

   return;

}

//------------------------------------------------------------------------------
template< >
int GetMethodData< MORTAR_WEIGHTS >( CouplingScheme const * cs )
{
   InterfacePairs const * const pairs = cs->getInterfacePairs();
   IndexType const numPairs = pairs->getNumPairs();

   ContactPlaneManager& cpManager = ContactPlaneManager::getInstance();
   parameters_t& parameters = parameters_t::getInstance();
   integer const dim = parameters.dimension;
   IndexType const numNodesPerFace = (dim == 3) ? 4 : 2;

   IndexType const masterId = cs->getMeshId1();
   IndexType const slaveId = cs->getMeshId2();

   ////////////////////////////////
   //                            //
   // compute single mortar gaps //
   //                            //
   ////////////////////////////////
   ComputeSingleMortarGaps( cs );

   int numRows = cs->getNumTotalNodes();
   static_cast<MortarData*>( cs->getMethodData() )->allocateMfemSparseMatrix( numRows );

   //////////////////////////////////////////////
   //                                          //
   // aggregate data to compute mortar weights //
   //                                          //
   //////////////////////////////////////////////

   // declare local variables to hold face nodal coordinates
   // and overlap vertex coordinates
   IndexType size = dim * numNodesPerFace;

   // projected coords
   real masterX_bar[ size ];
   real slaveX_bar[ size ];
   initRealArray( &masterX_bar[0], size, 0. );
   initRealArray( &slaveX_bar[0], size, 0. );

   real* overlapX; // [dim * cpManager.m_numPolyVert[cpID]];  

   int cpID = 0;
   for (IndexType kp = 0; kp < numPairs; ++kp)
   {
      InterfacePair pair = pairs->getInterfacePair(kp);

      if (!pair.inContact)
      {
         continue;
      }

      // get pair indices
      IndexType index1 = pair.pairIndex1;
      IndexType index2 = pair.pairIndex2;

      overlapX = new real[ dim * cpManager.m_numPolyVert[ cpID ]];
      initRealArray( overlapX, dim*cpManager.m_numPolyVert[cpID], 0. );

      // get projected face coordinates
      cpManager.getProjectedFaceCoords( cpID, 0, &masterX_bar[0] ); // face 0 = first face
      cpManager.getProjectedFaceCoords( cpID, 1, &slaveX_bar[0] ); // face 1 = second face

      // construct array of polygon overlap vertex coordinates
      cpManager.getContactPlaneOverlapVerts( cpID, cpManager.m_numPolyVert[cpID], &overlapX[0] );

      // instantiate contact surface element for purposes of computing 
      // mortar weights. Note, this uses projected face coords
      SurfaceContactElem elem( dim, &masterX_bar[0], &slaveX_bar[0], 
                               &overlapX[0],
                               numNodesPerFace, 
                               cpManager.m_numPolyVert[cpID],
                               masterId, slaveId, index1, index2 );

      // compute the mortar weights to be stored on the surface 
      // contact element struct. This must be done prior to computing nodal gaps
      elem.overlapArea = cpManager.m_area[ cpID ];

      ComputeMortarWeights( elem );

      elem.numActiveGaps = numNodesPerFace;

      // assemble mortar weight contributions sum_alpha int_alpha phi_a phi_b da.
      // Note: active slave nodes (i.e. active gaps) are checked in this routine.
      const EnforcementOptions& enforcement_options = const_cast<EnforcementOptions&>(cs->getEnforcementOptions());
      const SparseMode sparse_mode = enforcement_options.lm_implicit_options.sparse_mode;
      SLIC_ERROR_IF(sparse_mode == SparseMode::MFEM_ELEMENT_DENSE, 
         "GetMethodData<MORTAR_WEIGHTS>() MFEM_ELEMENT_DENSE " << 
         "Unassembled element dense matrix output not implemented.");
      static_cast<MortarData*>( cs->getMethodData() )->assembleMortarWts( elem, sparse_mode );

      ++cpID;

      delete [] overlapX;

   } // end loop over pairs to assemble mortar weights

   return 0;

} // end GetMethodData< MORTAR_WEIGHTS >()

//------------------------------------------------------------------------------

} // end namespace tribol
