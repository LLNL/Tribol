// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Tribol Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (MIT)

#ifndef SRC_MESH_MESHDATA_HPP_
#define SRC_MESH_MESHDATA_HPP_

#include "tribol/types.hpp"
#include "tribol/common/Parameters.hpp"
#include "tribol/utils/DataManager.hpp"

#include <ostream>

namespace tribol
{

/*!
 * \brief Struct to hold method specific nodal fields
 */
struct MeshNodalData 
{
   
   /// Default constructor
   MeshNodalData( );

   /// Destructor
   ~MeshNodalData() { }

   int m_numNodes;

   /////////////////////////
   // MORTAR NODAL FIELDS //
   /////////////////////////
   RealT * m_node_gap;            ///< scalar nodal gap (used on nonmortar mesh) 
   const RealT * m_node_pressure; ///< scalar nodal pressure (used on nonmortar mesh) 

   bool m_isGapComputed        {false}; ///< true if the nodal gaps have been computed
   bool m_is_node_gap_set      {false}; ///< true if nodal gap field is set
   bool m_is_node_pressure_set {false}; ///< true if nodal pressure field is set
   /////////////////////////

   bool m_is_velocity_set           {false}; ///< true if nodal velocities have been registered
   bool m_is_nodal_displacement_set {false}; ///< true if nodal displacements have been registered
   bool m_is_nodal_response_set     {false}; ///< true if the nodal responses have been registered

}; // end of struct MeshNodalData

/*!
 * \brief Struct to hold method specific element data
 */
struct MeshElemData 
{

   /// Default constructor
   MeshElemData();

   /// Destructor
   ~MeshElemData();

   int m_numCells;

   //////////////////////////////////////
   // PENALTY ENFORCEMENT ELEMENT DATA //
   //////////////////////////////////////
   const RealT * m_mat_mod;             ///< Bulk/Young's modulus for contact faces
   const RealT * m_thickness;           ///< Volume element thickness associated with each contact face

   RealT m_penalty_stiffness {0.};      ///< single scalar kinematic penalty stiffness for each mesh
   RealT m_penalty_scale {1.};          ///< scale factor applied to kinematic penalty only
   RealT m_rate_penalty_stiffness {0.}; ///< single scalar rate penalty stiffness for each mesh
   RealT m_rate_percent_stiffness {0.}; ///< rate penalty is percentage of gap penalty

   bool m_is_kinematic_constant_penalty_set {false}; ///< True if single kinematic constant penalty is set
   bool m_is_kinematic_element_penalty_set  {false}; ///< True if the element-wise kinematic penalty is set
   bool m_is_rate_constant_penalty_set      {false}; ///< True if the constant rate penalty is set
   bool m_is_rate_percent_penalty_set       {false}; ///< True if the rate percent penalty is set

   bool m_is_element_thickness_set          {false}; ///< True if element thickness is set

  /*!
   * \brief Checks if the kinematic penalty data is valid
   *
   * \param [in] pen_enfrc penalty enforcement option struct
   *
   * \Return True if the kinematic penalty option has valid data
   */
   bool isValidKinematicPenalty( PenaltyEnforcementOptions& pen_options );

  /*!
   * \brief Checks if the rate penalty data is valid
   *
   * \param [in] pen_enfrc to penalty enforcement option struct
   *
   * \Return True if the rate penalty option has valid data
   */
   bool isValidRatePenalty( PenaltyEnforcementOptions& pen_options ); ///< True if rate penalty option is valid

};

class MeshData
{
public:

 /*!
  * \brief Constructor 
  *
  */
  MeshData();

 /*!
  * \brief Destructor
  *
  */
  ~MeshData();

  InterfaceElementType m_elementType; ///< Type of interface element in mesh
  int m_lengthNodalData;              ///< Total number of elements in data arrays (used to create new arrays an index in using connectivity ids)
  int m_numSurfaceNodes;              ///< Total number of unique node ids in the surface connectivity (computed from MeshData::sortSurfaceNodeIds() )
  int m_numCells;                     ///< Total number of SURFACE cells in the mesh
  int m_numNodesPerCell;                 ///< Number of nodes per SURFACE cell based on cell type
  int m_dim;                          ///< Dimension of mesh
  int m_mesh_id;                       ///< Mesh Id associated with this data
  bool m_isValid;                     ///< True if the mesh is valid

  const RealT* m_positionX; ///< X-coordinates of nodes in mesh 
  const RealT* m_positionY; ///< Y-coordinates of nodes in mesh 
  const RealT* m_positionZ; ///< Z-coordinates of nodes in mesh 

  const RealT* m_dispX; ///< X-component of nodal displacements 
  const RealT* m_dispY; ///< Y-component of nodal displacements 
  const RealT* m_dispZ; ///< Z-component of nodal displacements 

  RealT * m_forceX; ///< X-component of nodal forces 
  RealT * m_forceY; ///< Y-component of nodal forces 
  RealT * m_forceZ; ///< Z-component of nodal forces 

  const RealT* m_velX; ///< X-component of nodal velocity 
  const RealT* m_velY; ///< Y-component of nodal velocity
  const RealT* m_velZ; ///< Z-component of nodal velocity

  RealT* m_nX; ///< X-component of outward unit face normals
  RealT* m_nY; ///< Y-component of outward unit face normals
  RealT* m_nZ; ///< Z-component of outward unit face normals

  RealT* m_cX; ///< X-component of vertex averaged cell centroids
  RealT* m_cY; ///< Y-component of vertex averaged cell centroids
  RealT* m_cZ; ///< Z-component of vertex averaged cell centroids

  RealT* m_faceRadius; ///< Face radius used in low level proximity check

  RealT* m_area; ///< Cell areas

  const IndexT* m_connectivity; ///< Cell connectivity arrays of length, m_numCells * m_numNodesPerCell

  IndexT* m_sortedSurfaceNodeIds; ///< List of sorted node ids from connectivity w/o duplicates, length = m_numSurfaceNodes.

  RealT* m_node_nX; ///< X-component of outward unit node normals
  RealT* m_node_nY; ///< Y-component of outward unit node normals
  RealT* m_node_nZ; ///< Z-component of outward unit node normals

  MeshNodalData m_nodalFields; ///< method specific nodal fields
  MeshElemData  m_elemData;    ///< method/enforcement specific element data 

public:

  /*!
   * \brief Checks for valid Lagrange multiplier enforcement data
   *
   */
   int checkLagrangeMultiplierData();

  /*!
   * \brief Checks for valid penalty enforcement data 
   *
   * \param [in] p_enfrc_options penalty enforcement options guiding check
   */
   int checkPenaltyData( PenaltyEnforcementOptions& p_enfrc_options );

  /// \brief Returns the number of cells in this Mesh instance
  int getNumberOfCells() const { return m_numCells; }


  /*!
   * \brief Computes the face normals and centroids for all faces in the mesh.
   *
   * \param [in] dim Dimension of the problem 
   * \return true if face calculations do not encounter errors or warnings
   * 
   * This routine accounts for warped faces by computing an average normal.
   */
  bool computeFaceData( int const dim );
  
  /*!
   * \brief Computes average nodal normals for use with mortar methods
   *
   * \note this routine computes average nodal normals for all nodes in the mesh.
   *
   * \param [in] dim Dimension of the problem
   */
  void computeNodalNormals( int const dim );

  /*!
   * \brief Get face node id in mesh 
   *
   * \param [in] faceId Id of face in mesh 
   * \param [in] localNodeId Local id (e.g. 0-3 for quad 4)
   *
   * \return Global mesh node id
   */
  int getFaceNodeId(int faceId, int localNodeId) const
     { return m_connectivity[m_numNodesPerCell*faceId+localNodeId]; }
  
  /*!
   *
   * \brief compute the approximate radius of the face's enclosing circle
   *
   * \param [in] faceId face id
   *
   * \return radius
   *
   * \note this routine finds the largest magnitude vector between a given 
   *  face vertex and the vertex averaged centroid and uses this as the 
   *  enclosing circle's radius. This is not necessarily the smallest 
   *  enclosing circle, but this computation is fast and serves the purpose 
   *  of the tribol proximity check.
   *
   */
  RealT computeFaceRadius( int faceId );
  
  /*!
   *
   * \brief compute the surface edge/segment length
   *
   * \param [in] edgeId edge id
   *
   * \return edge length
   *
   */
  RealT computeEdgeLength( int edgeId );
  
  /*!
   *
   * \brief returns pointer to array of stacked nodal coordinates for given face
   *
   * \param [in/out] coords pointer to an array of stacked (x,y,z) nodal coordinates
   *
   */
  void getFaceCoords( int const faceId, RealT * coords );

  /*!
   *
   * \brief returns pointer to array of stacked nodal velocities for given face
   *
   * \param [in/out] nodalVel pointer to an array of stacked (x,y,z) nodal velocities
   *
   */
  void getFaceNodalVelocities( int const faceId, RealT * nodalVel );

  /*!
   *
   * \brief returns pointer to array of stacked normal components 
   *
   * \param [in] faceId integer id of face
   * \param [in] dim dimension of the problem
   * \param [in/out] nrml pointer to array of stacked components of the face normal vector
   *
   */
  void getFaceNormal( int const faceId, int const dim, RealT * nrml );

  /// sorts unique surface node ids from connectivity and stores them on the mesh object in ascending order
  void sortSurfaceNodeIds();

  // Note, this routine is not used at the moment
  int localIdFromConn( const int connectivityId );
                    
  /*!
   * \brief Allocates storage space for face normals, centroids and areas.
   *
   * \note This routine assumes the m_numCells has been set to the number of 
   * mesh cells.
   * \param dimension The dimension of the mesh
   * \pre dimension is 2 or 3
   */
  void allocateArrays(int dimension);
  
  /*!
   * \brief Deallocates storage space for face normals, centroids and areas
   */
  void deallocateArrays();

  /// Prints information associated with this mesh to \a os
  void print(std::ostream& os) const;
};

using MeshManager = DataManager<MeshData>;

} // end namespace tribol

/// \a ostream operator to print a \a MeshData instance to \a os
std::ostream& operator<<(std::ostream& os, const tribol::MeshData& md);

#endif /* SRC_MESH_MESHDATA_HPP_ */
