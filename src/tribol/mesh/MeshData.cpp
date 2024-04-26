// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Tribol Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (MIT)

#include "tribol/mesh/MeshData.hpp"
#include "tribol/common/ExecModel.hpp"
#include "tribol/utils/Math.hpp"

#include <cmath> 
#include <iostream> 
#include <sstream>
#include <iomanip>
#include <fstream>

#include "axom/slic.hpp"
#include "axom/fmt.hpp"

namespace fmt = axom::fmt;

namespace tribol
{

//------------------------------------------------------------------------------
bool MeshElemData::isValidKinematicPenalty( PenaltyEnforcementOptions& pen_options )
{
  // Note, this routine is, and should be called only for non-null meshes
  KinematicPenaltyCalculation kin_calc = pen_options.kinematic_calculation;

  // check kinematic penalty calculation data
  if ( !in_range( kin_calc, NUM_KINEMATIC_PENALTY_CALCULATION) )
  {
    // warning already issued when penalty options were set. If penalty options 
    // were not set, the following else-if will catch this
    return false;
  }
  // a kinematic penalty should always be set. Right now Tribol does not support 
  // rate only enforcement
  else if ( !pen_options.is_kinematic_calculation_set() )
  {
    SLIC_WARNING( "MeshElemData::isValidKinematic(): kinematic penalty calculation data not set; " << 
                  "call tribol::setPenaltyOptions()." );
    return false;
  }

  switch (kin_calc)
  {
    case KINEMATIC_CONSTANT:
    {
      if ( !this->m_is_kinematic_constant_penalty_set )
      {
        SLIC_WARNING("MeshElemData::isValidKinematicPenalty(): " << 
                     "single stiffness penalty not set.");
        return false;
      }
      else if ( this->m_penalty_stiffness < pen_options.tiny_penalty )
      {
        SLIC_WARNING("MeshElemData::isValidKinematicPenalty(): " << 
                     "single stiffness penalty less than threshold (" << 
                     pen_options.tiny_penalty << "). Consider increasing " << 
                     "for your problem.");
        return false;
      }
      break;
    } // end case KINEMATIC_CONSTANT
    case KINEMATIC_ELEMENT:
    {
      if ( !this->m_is_kinematic_element_penalty_set )
      {
        SLIC_WARNING("MeshElemData::isValidKinematicPenalty(): " << 
                     "element-wise penalty data not set.");
        return false;
      }

      // check for positive material modulus and thickness values
      bool isValidMatMod = true;
      bool isValidElemThickness = true;
      for (int i=0; i<this->m_num_cells; ++i)
      {
        if (this->m_mat_mod[i] <= 0.)
        {
          isValidMatMod = false;
        }
        if (this->m_thickness[i] <= 0.)
        {
          isValidElemThickness = false;
        }

        SLIC_WARNING_IF(!isValidMatMod, "MeshElemData::isValidKinematicPenalty(): " << 
                        "invalid nonpositive element material modulus encountered.");

        SLIC_WARNING_IF(!isValidElemThickness, "MeshElemData::isValidKinematicPenalty(): " << 
                        "invalid nonpositive element thickness encountered.");

        if (!isValidMatMod || !isValidElemThickness)
        {
          return false;
        }
      } // end for loop over elements
      break;
    } // end case KINEMATIC_ELEMENT
    default:
      // no-op, quiet compiler
      break;
  } // end switch over kinematic calculation

  return true;
 
} // end MeshElemData::isValidKinematicPenalty()

//------------------------------------------------------------------------------
bool MeshElemData::isValidRatePenalty( PenaltyEnforcementOptions& pen_options )
{
  // Note, this method is and should only be called for non-null meshes
  RatePenaltyCalculation rate_calc = pen_options.rate_calculation; 

  // check rate penalty calculation data
  if ( !in_range( rate_calc, NUM_RATE_PENALTY_CALCULATION) )
  {
    // warning already issued when penalty options were set. If penalty options 
    // were not set, the following else-if will catch this
    return false;
  }
  // the rate_calc could be set to NONE and this boolean will be true
  else if ( !pen_options.is_rate_calculation_set() )
  {
    SLIC_WARNING( "MeshElemData::isValidRatePenalty(): rate penalty calculation data not set. " << 
                  "call tribol::setPenaltyOptions()." );
    return false;
  }

  switch (rate_calc)
  {
    case NO_RATE_PENALTY:
    {
      // double check that mesh booleans are consistent with rate penalty 
      // calculation option
      if (this->m_is_rate_constant_penalty_set)
      {
        this->m_is_rate_constant_penalty_set = false;
      }
      else if (this->m_is_rate_percent_penalty_set)
      {
        this->m_is_rate_percent_penalty_set = false;
      }
      break;
    } // end case NONE
    case RATE_CONSTANT:
    {
      if ( !this->m_is_rate_constant_penalty_set )
      {
        SLIC_WARNING("MeshElemData::isValidRatePenalty(): " << 
                     "constant rate penalty data not set.");
        return false;
      }
      else if ( this->m_rate_penalty_stiffness < pen_options.tiny_penalty )
      {
        SLIC_WARNING("MeshElemData::isValidRatePenalty(): " << 
                     "constant rate penalty less than threshold (" << 
                     pen_options.tiny_penalty << "). Consider increasing " << 
                     "for your problem.");
        return false;
      }
      break;
    } // end case RATE_CONSTANT
    case RATE_PERCENT:
    {
      if ( !this->m_is_rate_percent_penalty_set )
      {
        SLIC_WARNING("MeshElemData::isValidRatePenalty(): " << 
                     "percent rate penalty data not set.");
        return false;
      }
      else if ( this->m_rate_percent_stiffness < pen_options.tiny_penalty ||
                this->m_rate_percent_stiffness > (1.-pen_options.tiny_penalty) )
      {
        SLIC_WARNING("MeshElemData::isValidRatePenalty(): " << 
                     "rate percent penalty not in (0,1).");
        return false;
      }
      break;
    } // end case RATE_PERCENT
    default:
      // no-op, quiet compiler
      break;
  } // end switch on rate calculation

  return true;

} // end MeshElemData::isValidRatePenalty()

//------------------------------------------------------------------------------

/////////////////////////////////
//                             //
// Routines for MeshData class //
//                             //
/////////////////////////////////
MeshData::MeshData( IndexT mesh_id, InterfaceElementType element_type,
    IndexT num_nodes, IndexT num_elements )
  : m_mesh_id( mesh_id )
  , m_element_type( element_type )
  , m_dim( getDimFromElementType() )
  , m_num_nodes( num_nodes )
  , m_num_elements( num_elements )
  , m_is_valid( true )
{}


template <MemorySpace MSPACE>
MeshDataBySpace<MSPACE>::MeshDataBySpace( IndexT mesh_id, IndexT num_elements, IndexT num_nodes,
                    const IndexT* connectivity, InterfaceElementType element_type,
                    const RealT* x, const RealT* y, const RealT* z )
  : MeshData( mesh_id, element_type, num_nodes, num_elements )
  , m_position( createNodalVector(x, y, z) )
  , m_connectivity( createConnectivity(num_elements, connectivity) )
{
  // mesh verification
  if (num_elements > 0)
  {
    if (m_dim == 2 && (x == nullptr || y == nullptr))
    {
      SLIC_WARNING_ROOT("tribol::MeshData(): pointer to x and/or y-component " << 
                        "mesh coordinate array is a null pointer " <<
                        "for mesh id " << m_mesh_id << ".");
      m_is_valid = false;
    }
    else if (m_dim == 3 && (x == nullptr || y == nullptr || z == nullptr))
    {
      SLIC_WARNING_ROOT("tribol::MeshData(): pointer to x, y, and/or z-component " << 
                        "mesh coordinate array is a null pointer " <<
                        "for mesh id " << m_mesh_id << ".");
      m_is_valid = false;
    }
    if (connectivity == nullptr)
    {
      SLIC_WARNING_ROOT("tribol::MeshData(): pointer to mesh connectivity is " <<
                        "a null pointer for mesh id " << m_mesh_id << ".");
    }
  }

  // Find unique surface node ids
  if (MSPACE == MemorySpace::Host)
  {
    if (num_elements > 0)
    {
      sortSurfaceNodeIds();
    }
  }

  getElementData().m_num_cells = num_elements;
  getNodalFields().m_num_nodes = num_nodes;
}

//------------------------------------------------------------------------------
template <MemorySpace MSPACE>
void MeshDataBySpace<MSPACE>::setPosition( const RealT* x,
                                           const RealT* y,
                                           const RealT* z )
{
  m_position = createNodalVector(x, y, z);
}

//------------------------------------------------------------------------------
template <MemorySpace MSPACE>
void MeshDataBySpace<MSPACE>::setDisplacement( const RealT* ux,
                                               const RealT* uy,
                                               const RealT* uz )
{
  m_disp = createNodalVector(ux, uy, uz);
}

//------------------------------------------------------------------------------
template <MemorySpace MSPACE>
void MeshDataBySpace<MSPACE>::setVelocity( const RealT* vx,
                                           const RealT* vy,
                                           const RealT* vz )
{
  m_vel = createNodalVector(vx, vy, vz);
}

//------------------------------------------------------------------------------
template <MemorySpace MSPACE>
void MeshDataBySpace<MSPACE>::setResponse( RealT* rx,
                                           RealT* ry,
                                           RealT* rz )
{
  m_response = createNodalVector(rx, ry, rz);
}

//------------------------------------------------------------------------------
int MeshData::getDimFromElementType()
{
  switch (m_element_type)
  {
    case LINEAR_EDGE:
    {
      return 2;
    }
    case LINEAR_TRIANGLE:
    case LINEAR_QUAD:
    {
      return 3;
    }
    default:
    {
      SLIC_ERROR_ROOT("Unsupported element type for a contact mesh.");
      return 0;
    }
  }
}

//------------------------------------------------------------------------------
template <MemorySpace MSPACE>
ArrayViewT<const IndexT, 2, MSPACE> MeshDataBySpace<MSPACE>::createConnectivity(IndexT num_elements, const IndexT* connectivity)
{
  switch (m_element_type)
  {
    case LINEAR_EDGE:
    {
      return ArrayViewT<const IndexT, 2, MSPACE>(connectivity, {num_elements, 2});
    }
    case LINEAR_TRIANGLE:
    {
      return ArrayViewT<const IndexT, 2, MSPACE>(connectivity, {num_elements, 3});
    }
    case LINEAR_QUAD:
    {
      return ArrayViewT<const IndexT, 2, MSPACE>(connectivity, {num_elements, 4});
    }
    default:
    {
      SLIC_ERROR_ROOT("Unsupported element type for a contact mesh.");
      return ArrayViewT<const IndexT, 2, MSPACE>(connectivity, {num_elements, 0});
    }
  }
}

//------------------------------------------------------------------------------
template <MemorySpace MSPACE>
void MeshDataBySpace<MSPACE>::sortSurfaceNodeIds()
{
  ArrayT<IndexT> sorted_conn(0, m_connectivity.size());
  for (auto node_id : m_connectivity)
  {
    sorted_conn.push_back(node_id);
  }

  bubble_sort(sorted_conn.data(), sorted_conn.size());

  // count number of duplicate entries
  int num_dup = 0;
  for (IndexT i{1}; i < sorted_conn.size(); ++i)
  {
    if (sorted_conn[i] == sorted_conn[i-1])
    {
      ++num_dup;
    }
  }

  // compute number of unique integer ids
  int unique_size = sorted_conn.size() - num_dup;

  SLIC_ERROR_IF(unique_size <= 0, "MeshData::sortSurfaceNodeIds(): " << 
                "invalid connectivity array; " <<
                "only single unique id in connectivity array.");

  // allocate array to store unique, sorted node ids on mesh object
  m_sorted_surface_node_ids = ArrayT<IndexT>(0, unique_size);

  // populate sorted node id list
  m_sorted_surface_node_ids.push_back(sorted_conn[0]);
  for (IndexT i{1}; i < sorted_conn.size(); ++i)
  {
    if ( sorted_conn[i] != sorted_conn[i-1] )
    { 
      m_sorted_surface_node_ids.push_back(sorted_conn[i]);
    }
  }

  return;
} // end MeshData::sortSurfaceNodeIds()

//------------------------------------------------------------------------------
template <MemorySpace MSPACE>
bool MeshDataBySpace<MSPACE>::computeFaceData()
{
  bool faceDataOk = true;
  constexpr RealT nrmlMagTol = 1.0e-15;
  const IndexT num_elements = numberOfElements();

  // allocate m_c
  ArrayT<ArrayT<RealT, 1, MSPACE>, 1, MemorySpace::Dynamic> c_host(m_dim, m_dim);
  for (auto& c_dim : c_host)
  {
    c_dim = ArrayT<RealT, 1, MSPACE>(num_elements, num_elements);
  }
  m_c = std::move(c_host);

  // allocate face_radius
  m_face_radius = ScalarArray<RealT, MSPACE>(num_elements, num_elements);

  // allocate m_n
  ArrayT<ArrayT<RealT, 1, MSPACE>, 1, MemorySpace::Dynamic> n_host(m_dim, m_dim);
  for (auto& n_dim : n_host)
  {
    n_dim = ArrayT<RealT, 1, MSPACE>(num_elements, num_elements);
  }
  m_n = std::move(n_host);

  // allocate m_area
  m_area = ScalarArray<RealT, MSPACE>(num_elements, num_elements);

  // loop over all elements in the mesh
  auto num_nodes_per_elem = numberOfNodesPerElement();
  auto c_data = m_c.data();
  auto position_data = m_position.data();
  auto dim = m_dim;
  forAllExec<toExecutionMode<MSPACE>::value>(num_elements, 
    [dim, num_nodes_per_elem, c_data, position_data] TRIBOL_HOST_DEVICE (const IndexT i) {

     // compute the vertex average centroid. This will lie in the 
     // plane of the face for planar faces, and will be used as 
     // an approximate centroid for warped faces, both in 3D.

     // loop over the nodes per element
     for (int j=0; j<num_nodes_per_elem; ++j) {
        auto nodeId = 1;//getGlobalNodeId(i, j);
        // always compute the x and y components for 2D and 3D
        c_data[0][i] += position_data[0][ nodeId ];
        c_data[1][i] += position_data[1][ nodeId ];
     } // end loop over nodes
     
     RealT fac = 1.0 / num_nodes_per_elem;
     c_data[0][i] = fac * c_data[0][i];
     c_data[1][i] = fac * c_data[1][i];

     if (dim == 3) {
        for (int j=0; j<num_nodes_per_elem; ++j) {
           auto nodeId = 1;//getGlobalNodeId(i, j);
           c_data[2][i] += position_data[2][ nodeId ];
        } // end loop over nodes
        c_data[2][i] = fac * c_data[2][i];
     } // end if-dim

    //  // compute the outward facing normal
    //  if (m_dim == 2) {

    //     // the 2D calculation over a 2-node, 1D segment assumes a 
    //     // counter-clockwise ordering of the quad4 area element
    //     // to which the 1D line segment belongs. This is to properly 
    //     // orient the normal outward
    //     auto nodeId = getGlobalNodeId(i, 0);
    //     auto nextNodeId = getGlobalNodeId(i, 1);
    //     RealT lambdaX = position_data[0][ nextNodeId ] - position_data[0][ nodeId ];
    //     RealT lambdaY = position_data[1][ nextNodeId ] - position_data[1][ nodeId ];
   
    //     m_n[0][i] = lambdaY;
    //     m_n[1][i] = -lambdaX;

    //     // compute the length of the segment
    //     m_area[i] = magnitude( lambdaX, lambdaY );

    //     // normalize normal vector
    //     auto mag = magnitude( m_n[0][i], m_n[1][i] );
    //     auto invMag = nrmlMagTol;
    //     if (mag >= nrmlMagTol) {
    //        invMag = 1.0 / mag;
    //     } else {
    //        //faceDataOk = false;
    //     }
    //     m_n[0][i] *= invMag;
    //     m_n[1][i] *= invMag;

    //  }
        
    //  // compute face radius for both 2D and 3D. For 2D, this is a duplicate of 
    //  // the length, but allows for some uniformity in accessing this value 
    //  // for tolerance ratios
    //  m_face_radius[i] = this->computeFaceRadius(i);

    //  if (m_dim == 3) {

    //     // this method of computing an outward unit normal breaks the 
    //     // face into triangular pallets by connecting two consecutive 
    //     // nodes with the approximate centroid.
    //     // The average outward unit normal for the face is the average of 
    //     // those of the pallets. This is exact for non-warped faces. To 
    //     // compute the pallet normal, you only need edge vectors for the 
    //     // pallet. These are constructed from the face centroid and the face 
    //     // edge's first node and the face edge's two nodes

    //     // declare triangle edge vector components and normal components
    //     RealT vX1, vY1, vZ1;
    //     RealT vX2, vY2, vZ2;
    //     RealT nX, nY, nZ; 

    //     // loop over num_nodes_per_elem-1 element edges and compute pallet normal
    //     for (int j=0; j<(num_nodes_per_elem-1); ++j) 
    //     {
    //        auto nodeId = getGlobalNodeId(i, j);
    //        auto nextNodeId = getGlobalNodeId(i, j+1);
    //        // first triangle edge vector between the face's two 
    //        // edge nodes
    //        vX1 = position_data[0][ nextNodeId ] - position_data[0][ nodeId ];
    //        vY1 = position_data[1][ nextNodeId ] - position_data[1][ nodeId ];
    //        vZ1 = position_data[2][ nextNodeId ] - position_data[2][ nodeId ];
          
    //        // second triangle edge vector between the face centroid 
    //        // and the face edge's first node
    //        vX2 = m_c[0][i] - position_data[0][ nodeId ];
    //        vY2 = m_c[1][i] - position_data[1][ nodeId ];
    //        vZ2 = m_c[2][i] - position_data[2][ nodeId ];

    //        // compute the contribution to the pallet normal as v1 x v2. Sum 
    //        // these into the face normal component variables stored on the mesh data 
    //        // object
    //        nX = (vY1 * vZ2) - (vZ1 * vY2);
    //        nY = (vZ1 * vX2) - (vX1 * vZ2);
    //        nZ = (vX1 * vY2) - (vY1 * vX2);

    //        // sum the normal component contributions into the component variables
    //        m_n[0][i] += nX;
    //        m_n[1][i] += nY;
    //        m_n[2][i] += nZ;

    //        // half the magnitude of the computed normal is the pallet area. Note: this is exact 
    //        // for planar faces and approximate for warped faces. Face areas are used in a general 
    //        // sense to create a face-overlap tolerance
    //        m_area[i] += 0.5 * magnitude( nX, nY, nZ );
    //     }

    //     // compute the pallet normal contribution for the last pallet
    //     auto nodeId = getGlobalNodeId(i, 0);
    //     auto nextNodeId = m_connectivity(i, num_nodes_per_elem - 1);
    //     vX1 = position_data[0][ nodeId ] - position_data[0][ nextNodeId ];
    //     vY1 = position_data[1][ nodeId ] - position_data[1][ nextNodeId ];
    //     vZ1 = position_data[2][ nodeId ] - position_data[2][ nextNodeId ];

    //     vX2 = m_c[0][i] - position_data[0][ nextNodeId ];
    //     vY2 = m_c[1][i] - position_data[1][ nextNodeId ];
    //     vZ2 = m_c[2][i] - position_data[2][ nextNodeId ];
        
    //     nX = (vY1 * vZ2) - (vZ1 * vY2);
    //     nY = (vZ1 * vX2) - (vX1 * vZ2);
    //     nZ = (vX1 * vY2) - (vY1 * vX2);

    //     // sum the normal component contributions into the component variables
    //     m_n[0][i] += nX;
    //     m_n[1][i] += nY;
    //     m_n[2][i] += nZ;

    //     // half the magnitude of the computed normal is the pallet area. Note: this is exact 
    //     // for planar faces and approximate for warped faces. Face areas are used in a general 
    //     // sense to create a face-overlap tolerance
    //     m_area[i] += 0.5 * magnitude( nX, nY, nZ );

    //     // multiply the pallet normal components by fac to obtain avg.
    //     m_n[0][i] = fac * m_n[0][i];
    //     m_n[1][i] = fac * m_n[1][i];
    //     m_n[2][i] = fac * m_n[2][i];

    //     // compute the magnitude of the average pallet normal
    //     auto mag = magnitude(m_n[0][i], m_n[1][i], m_n[2][i] );
    //     auto invMag = nrmlMagTol;
    //     if (mag >= nrmlMagTol) {
    //        invMag = 1.0 / mag;
    //     } else {
    //        //faceDataOk = false;
    //     }

    //     // normalize the average normal
    //     m_n[0][i] *= invMag;
    //     m_n[1][i] *= invMag;
    //     m_n[2][i] *= invMag;

    //  } // end if (dim == 3)

  }); // end element loop

  SLIC_WARNING_IF(!faceDataOk, 
      fmt::format("There are faces with a normal magnitude less than tolerance ({:e}).", nrmlMagTol));

  return faceDataOk; 

} // end MeshData::computeFaceData()

//------------------------------------------------------------------------------
template <MemorySpace MSPACE>
TRIBOL_HOST_DEVICE RealT MeshDataBySpace<MSPACE>::computeFaceRadius( int faceId ) 
{
   // loop over nodes of the face and determine the maximum 
   // "link" vector from the ith node to the face center
   RealT sqrRadius = 0.0;
   for (int i=0; i<numberOfNodesPerElement(); ++i) {
      const int nodeId = getGlobalNodeId(faceId, i);
      RealT lvx = m_position[0][nodeId] - m_c[0][faceId];
      RealT lvy = m_position[1][nodeId] - m_c[1][faceId];

      RealT lvz;
      if (m_dim == 3) { // for 3D
         lvz = m_position[2][nodeId] - m_c[2][faceId];
      }
      else {
         lvz = 0.0;
      }

      RealT sqrLinkMag = lvx * lvx + lvy * lvy + lvz * lvz;   
     
      if (sqrLinkMag > sqrRadius) {
         sqrRadius = sqrLinkMag;
      }
   } 

   return sqrt(sqrRadius);

} // end MeshData::computeFaceRadius()

//------------------------------------------------------------------------------
template <MemorySpace MSPACE>
RealT MeshDataBySpace<MSPACE>::computeEdgeLength( int faceId ) 
{
   // compute the length of the edge as the magnitude of 
   // the vector defined between the two edge vertices
   const int nodeId1 = getGlobalNodeId( faceId, 0 );
   const int nodeId2 = getGlobalNodeId( faceId, 1 );
   RealT lvx = m_position[0][nodeId2] - m_position[0][nodeId1];
   RealT lvy = m_position[1][nodeId2] - m_position[1][nodeId1];

   RealT len = magnitude( lvx, lvy );

   return len;

} // end MeshData::computeEdgeLength()

//------------------------------------------------------------------------------
template <MemorySpace MSPACE>
void MeshDataBySpace<MSPACE>::getFaceCoords( int const faceId, RealT * coords )
{

   for (IndexT a=0; a<numberOfNodesPerElement(); ++a)
   {
     IndexT nodeId = getGlobalNodeId(faceId, a);

     coords[m_dim * a]     = m_position[0][ nodeId ];
     coords[m_dim * a + 1] = m_position[1][ nodeId ];

     if (m_dim == 3)
     {
        coords[m_dim * a + 2] = m_position[2][ nodeId ];
     }
   }

   return;

}  // end MeshData::getFaceCoords()

//------------------------------------------------------------------------------
template <MemorySpace MSPACE>
void MeshDataBySpace<MSPACE>::getFaceNodalVelocities( int const faceId, RealT * nodalVel )
{
   if (m_vel.empty())
   {
      SLIC_ERROR("MeshData::getFaceNodalVelocities(): velocity arrays on each " <<
                 "mesh must be registered with Tribol.");
   }

   for (IndexT a=0; a<numberOfNodesPerElement(); ++a)
   {
     IndexT nodeId = getGlobalNodeId(faceId, a);

     nodalVel[m_dim * a]     = m_vel[0][ nodeId ];
     nodalVel[m_dim * a + 1] = m_vel[1][ nodeId ];

     if (m_dim == 3)
     {
        nodalVel[m_dim * a + 2] = m_vel[2][ nodeId ];
     }
   }

   return;

}  // end MeshData::getFaceNodalVelocities()

//------------------------------------------------------------------------------
template <MemorySpace MSPACE>
void MeshDataBySpace<MSPACE>::computeNodalNormals( int const dim )
{
   int * numFaceNrmlsToNodes;
   if (this->numberOfElements() > 0)
   {
      // check to make sure face normals have been computed with 
      // a call to computeFaceData
      if (m_n.empty())
      {
         SLIC_ERROR("MeshData::computeNodalNormals: required face normals not computed.");
      }

      m_node_n = VectorArray<RealT, MSPACE>(m_dim, m_dim);
      for (IndexT i{0}; i < m_dim; ++i)
      {
        m_node_n[i] = ArrayT<RealT, 1, MSPACE>(numberOfNodes(), numberOfNodes());
      }

      // allocate space for nodal normal array
      int size = numberOfNodes(); 
      // allocate scratch array to hold number of faces whose 
      // normals contribute to a given node. Most of the time 
      // this will be four face normals contributing to an 
      // averaged nodal normal for linear quad elements, but 
      // we want to handle arbitrary meshes and edge cases
      allocIntArray( &numFaceNrmlsToNodes, size, 0 );
   } // end if-check on null mesh

   // loop over elements
   for (int i=0; i<this->numberOfElements(); ++i)
   {
      // loop over element nodes
      for (int j=0; j<this->numberOfNodesPerElement(); ++j)
      {

         // SRW: note the connectivity array must be local to the mesh for indexing into 
         // the mesh nodal normal array. If it is not, then nodeId will access some other
         // piece of memory and there may be a memory issue when numFaceNrmlsToNodes is deleted
         // at the end of this routine.
         int nodeId = getGlobalNodeId(i, j);
         m_node_n[0][ nodeId ] += this->m_n[0][ i ]; // m_n[0][i] is the ith face normal x-component
         m_node_n[1][ nodeId ] += this->m_n[1][ i ]; // see above...

         // increment face-normal-to-node contribution counter
         ++numFaceNrmlsToNodes[ nodeId ];

      } // end loop over element nodes

      // populate z-component for 3D problems
      if (dim == 3)
      {

         // repeat loop over element nodes for z-component
         for (int k=0; k<this->numberOfNodesPerElement(); ++k)
         {
            int nodeId = getGlobalNodeId(i, k); 
            m_node_n[2][ nodeId ] += this->m_n[2][ i ];
         } // end loop over element nodes

      } // end if (dim == 3)

   } // end loop over elements

   // average the nodal normals
   if (this->numberOfElements() > 0)
   {
      for (int i=0; i<numberOfNodes(); ++i)
      {
         m_node_n[0][i] /= numFaceNrmlsToNodes[i];
         m_node_n[1][i] /= numFaceNrmlsToNodes[i];
         if (dim == 3)
         {
            m_node_n[2][i] /= numFaceNrmlsToNodes[i];
         }
      } // end loop over nodes
   }

   // normalize the nodal normals
   if (this->numberOfElements() > 0)
   {
      if (dim == 3)
      {
         for (int i=0; i<numberOfNodes(); ++i)
         {
            RealT mag = magnitude( m_node_n[0][i], m_node_n[1][i], m_node_n[2][i] );
            m_node_n[0][ i ] /= mag;
            m_node_n[1][ i ] /= mag;
            m_node_n[2][ i ] /= mag;
         }
      }
      else 
      {
         for (int i=0; i<numberOfNodes(); ++i)
         {
            RealT mag = magnitude( m_node_n[0][i], m_node_n[1][i] );
            m_node_n[0][ i ] /= mag;
            m_node_n[1][ i ] /= mag;
         }
      }

      delete [] numFaceNrmlsToNodes;
   } // end if-check on null mesh

   return;
} // end MeshData::computeNodalNormals()

//------------------------------------------------------------------------------
template <MemorySpace MSPACE>
void MeshDataBySpace<MSPACE>::getFaceNormal( int const faceId, int const dim, RealT * nrml )
{
   nrml[0] = this->m_n[0][ faceId ];
   nrml[1] = this->m_n[1][ faceId ];
   
   if (dim == 3)
   {
      nrml[2] = this->m_n[2][ faceId ];
   }

   return;

} // end MeshData::getFaceNormal()

//------------------------------------------------------------------------------
int MeshData::checkLagrangeMultiplierData()
{
   int err = 0;
   if (this->numberOfElements()>0)
   {
      if (!m_nodal_fields.m_is_node_gap_set ||
          !m_nodal_fields.m_is_node_pressure_set)
      {
         err = 1;
      }
   } // end if-non-null mesh
   return err; 
}
//------------------------------------------------------------------------------
int MeshData::checkPenaltyData( PenaltyEnforcementOptions& p_enfrc_options )
{
   int err = 0;
   if (this->numberOfElements()>0)
   {
      PenaltyConstraintType constraint_type = p_enfrc_options.constraint_type;
      // switch over penalty enforcement options and check for required data
      switch (constraint_type)
      {
         case KINEMATIC:
         {
            if (!m_element_data.isValidKinematicPenalty( p_enfrc_options ))
            {
               err = 1;
            }
            break;
         } // end KINEMATIC case

         case KINEMATIC_AND_RATE:
         {
            if (!m_element_data.isValidKinematicPenalty( p_enfrc_options ))
            {
               err = 1;
            }
            if (!m_element_data.isValidRatePenalty( p_enfrc_options ))
            {
               err = 1;
            }
            if (!m_nodal_fields.m_is_velocity_set)
            {
               SLIC_WARNING("Nodal velocities not set or null pointers; please set for " << 
                            "use with gap rate penalty enforcement.");
               err = 1;
            }
            break;
         } // end case KINEMATIC_AND_RATE
         default:
            // no-op, quiet compiler
            break;
      } // end switch over constraint types
   } // end if-non-null mesh

   return err;
} // end MeshData::checkPenaltyData()

//------------------------------------------------------------------------------
template <MemorySpace MSPACE>
void MeshDataBySpace<MSPACE>::print(std::ostream& os) const
{
   const int num_verts = numberOfNodes();
   const int num_elem = numberOfElements();

   if(num_verts <= 0)
   {
      os << "{}";
      return;
   }

   os << "{\n";
   os << fmt::format("  verts ({}) {{",num_verts);
   // positions
   os << fmt::format("\n\tx: {}", fmt::join(m_position[0].data(), m_position[0].data()+num_verts, ", "));
   os << fmt::format("\n\ty: {}", fmt::join(m_position[1].data(), m_position[1].data()+num_verts, ", "));
   if(m_dim == 3)
   {  
      os << fmt::format("\n\tz: {}", fmt::join(m_position[2].data(), m_position[2].data()+num_verts, ", "));
   }
   // contact response (force)
   if( !m_response.empty() )
   {
      os << fmt::format("\n\tfx: {}", fmt::join(m_response[0].data(), m_response[0].data()+num_verts, ", "));
      os << fmt::format("\n\tfy: {}", fmt::join(m_response[1].data(), m_response[1].data()+num_verts, ", "));
      if(m_dim == 3)
      {  
         os << fmt::format("\n\tfz: {}", fmt::join(m_response[2].data(), m_response[2].data()+num_verts, ", "));
      }
   }
   os << "\n  }";

   os << fmt::format("\n  elems ({}) {{", num_elem);  
   
   if( !m_connectivity.empty() )
   {
      os << fmt::format("\n\tconnectivity: {{ {} }}", fmt::join(m_connectivity.data(), m_connectivity.data()+(num_elem*numberOfNodesPerElement()), ", "));
   }

   // normals
   if( !m_n.empty() )
   {
      os << fmt::format("\n\tnx: {}", fmt::join(m_n[0].data(), m_n[0].data()+num_elem, ", "));
      os << fmt::format("\n\tny: {}", fmt::join(m_n[1].data(), m_n[1].data()+num_elem, ", "));
      if(m_dim == 3)
      {  
         os << fmt::format("\n\tnz: {}", fmt::join(m_n[2].data(), m_n[2].data()+num_elem, ", "));
      }
   }
   os << "\n  }";

   os << "\n}";
}

//------------------------------------------------------------------------------
template class MeshDataBySpace<MemorySpace::Dynamic>;
#ifdef TRIBOL_USE_UMPIRE
template class MeshDataBySpace<MemorySpace::Host>;
template class MeshDataBySpace<MemorySpace::Device>;
template class MeshDataBySpace<MemorySpace::Unified>;
#endif

} // end tribol namespace

std::ostream& operator<<(std::ostream& os, const tribol::MeshData& md)
{
   md.print(os);
   return os;
}
