// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Tribol Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (MIT)

#ifndef SRC_MESH_MESHDATA_HPP_
#define SRC_MESH_MESHDATA_HPP_

#include <ostream>

#include "tribol/common/ArrayTypes.hpp"
#include "tribol/common/LoopExec.hpp"
#include "tribol/common/Parameters.hpp"
#include "tribol/utils/DataManager.hpp"

namespace tribol
{

/*!
 * \brief Struct to hold method specific nodal fields
 */
struct MeshNodalData 
{
  int m_num_nodes;

  /////////////////////////
  // MORTAR NODAL FIELDS //
  /////////////////////////
  ArrayViewT<RealT>       m_node_gap;      ///< scalar nodal gap (used on nonmortar mesh) 
  ArrayViewT<const RealT> m_node_pressure; ///< scalar nodal pressure (used on nonmortar mesh) 

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
  int m_num_cells;

  //////////////////////////////////////
  // PENALTY ENFORCEMENT ELEMENT DATA //
  //////////////////////////////////////
  ArrayViewT<const RealT> m_mat_mod;   ///< Bulk/Young's modulus for contact faces
  ArrayViewT<const RealT> m_thickness; ///< Volume element thickness associated with each contact face

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
  class Viewer
  {
  public:
    Viewer(MeshData& mesh);

    TRIBOL_HOST_DEVICE IndexT meshId() const { return m_mesh_id; }
    TRIBOL_HOST_DEVICE InterfaceElementType getElementType() const { return m_element_type; }
    TRIBOL_HOST_DEVICE MemorySpace getMemorySpace() const { return m_mem_space; }
    TRIBOL_HOST_DEVICE int getAllocatorId() const { return m_allocator_id; }

    TRIBOL_HOST_DEVICE MeshNodalData& getNodalFields() { return m_nodal_fields; }
    TRIBOL_HOST_DEVICE const MeshNodalData& getNodalFields() const { return m_nodal_fields; }
    TRIBOL_HOST_DEVICE MeshElemData& getElementData() { return m_element_data; }
    TRIBOL_HOST_DEVICE const MeshElemData& getElementData() const { return m_element_data; }

    TRIBOL_HOST_DEVICE int spatialDimension() const { return m_position.size(); }
    TRIBOL_HOST_DEVICE IndexT numberOfNodes() const { return m_num_nodes; }
    TRIBOL_HOST_DEVICE IndexT numberOfElements() const { return m_connectivity.shape()[0]; }
    TRIBOL_HOST_DEVICE IndexT numberOfNodesPerElement() const { return m_connectivity.shape()[1]; }
    TRIBOL_HOST_DEVICE IndexT getGlobalNodeId(IndexT element_id, IndexT local_node_id) const
    {
      return m_connectivity(element_id, local_node_id);
    }

    TRIBOL_HOST_DEVICE const MultiViewArrayView<const RealT>& getPosition() const
    {
      return m_position;
    }
    
    TRIBOL_HOST_DEVICE bool hasDisplacement() const { return !m_disp.empty(); }
    TRIBOL_HOST_DEVICE const MultiViewArrayView<const RealT>& getDisplacement() const
    {
      return m_disp;
    }

    TRIBOL_HOST_DEVICE bool hasVelocity() const { return !m_vel.empty(); }
    TRIBOL_HOST_DEVICE const MultiViewArrayView<const RealT>& getVelocity() const
    {
      return m_vel;
    }
    
    TRIBOL_HOST_DEVICE bool hasResponse() const { return !m_response.empty(); }
    TRIBOL_HOST_DEVICE const MultiViewArrayView<RealT>& getResponse() const
    {
      return m_response;
    }
    
    TRIBOL_HOST_DEVICE bool hasNodalNormals() const { return !m_node_n.empty(); }
    TRIBOL_HOST_DEVICE const Array2DView<RealT>& getNodalNormals() const
    {
      return m_node_n;
    }

    TRIBOL_HOST_DEVICE bool hasElementCentroids() const { return !m_c.empty(); }
    TRIBOL_HOST_DEVICE const Array2DView<RealT>& getElementCentroids() const
    {
      return m_c;
    }

    TRIBOL_HOST_DEVICE bool hasElementNormals() const { return !m_n.empty(); }
    TRIBOL_HOST_DEVICE const Array2DView<RealT>& getElementNormals() const
    {
      return m_n;
    }

    TRIBOL_HOST_DEVICE bool hasFaceRadii() const { return !m_face_radius.empty(); }
    TRIBOL_HOST_DEVICE const Array1DView<RealT>& getFaceRadii() const
    {
      return m_face_radius;
    }

    TRIBOL_HOST_DEVICE bool hasElementAreas() const { return !m_area.empty(); }
    TRIBOL_HOST_DEVICE const Array1DView<RealT>& getElementAreas() const
    {
      return m_area;
    }

    TRIBOL_HOST_DEVICE const Array2DView<const IndexT>& getConnectivity() const 
    {
      return m_connectivity;
    }
  
    /*!
    *
    * \brief returns pointer to array of stacked nodal coordinates for given face
    *
    * \param [in] face_id integer id of face
    * \param [in/out] coords pointer to an array of stacked (x,y,z) nodal coordinates
    *
    */
    TRIBOL_HOST_DEVICE void getFaceCoords( IndexT face_id, RealT* coords ) const;

    /*!
    *
    * \brief returns pointer to array of stacked nodal velocities for given face
    *
    * \param [in] face_id integer id of face
    * \param [in/out] nodalVel pointer to an array of stacked (x,y,z) nodal velocities
    *
    */
    TRIBOL_HOST_DEVICE void getFaceVelocities( IndexT face_id, RealT* vels ) const;

    /*!
    *
    * \brief returns pointer to array of stacked normal components 
    *
    * \param [in] face_id integer id of face
    * \param [in/out] nrml pointer to array of stacked components of the face normal vector
    *
    */
    TRIBOL_HOST_DEVICE void getFaceNormal( IndexT face_id, RealT* nrml ) const;
    
  private:
    const IndexT m_mesh_id;
    const InterfaceElementType m_element_type;
    const IndexT m_num_nodes;

    const MemorySpace m_mem_space;
    const int m_allocator_id;

    const MultiViewArrayView<const RealT> m_position;
    const MultiViewArrayView<const RealT> m_disp;
    const MultiViewArrayView<const RealT> m_vel;
    const MultiViewArrayView<RealT> m_response;

    const Array2DView<RealT> m_node_n;

    const Array2DView<const IndexT> m_connectivity;

    const Array2DView<RealT> m_c;
    const Array2DView<RealT> m_n;
    const ArrayViewT<RealT> m_face_radius;
    const ArrayViewT<RealT> m_area;
    
    MeshNodalData m_nodal_fields; ///< method specific nodal fields
    MeshElemData  m_element_data; ///< method/enforcement specific element data
  };

  /*!
  * \brief Constructor 
  *
  */
  MeshData( IndexT mesh_id, IndexT num_elements, IndexT num_nodes,
                   const IndexT* connectivity, InterfaceElementType element_type,
                   const RealT* x, const RealT* y, const RealT* z,
                   MemorySpace mem_space );

  InterfaceElementType getElementType() const { return m_element_type; }
  int spatialDimension() const { return m_dim; }
  MemorySpace getMemorySpace() const { return m_mem_space; }
  int getAllocatorId() const { return m_allocator_id; }
  void updateAllocatorId(int allocator_id ) { m_allocator_id = allocator_id; }

  bool& isMeshValid() { return m_is_valid; }

  MeshNodalData& getNodalFields() { return m_nodal_fields; }
  MeshElemData& getElementData() { return m_element_data; }
  const MeshElemData& getElementData() const { return m_element_data; }
  
  IndexT numberOfNodes() const { return m_num_nodes; }
  IndexT numberOfElements() const { return m_connectivity.shape()[0]; }
  IndexT numberOfNodesPerElement() const { return m_connectivity.shape()[1]; }
  IndexT getGlobalNodeId(IndexT element_id, IndexT local_node_id) const
  {
    return m_connectivity(element_id, local_node_id);
  }

  void setPosition( const RealT* x,
                    const RealT* y,
                    const RealT* z );

  void setDisplacement( const RealT* ux,
                        const RealT* uy,
                        const RealT* uz );

  void setVelocity( const RealT* vx,
                    const RealT* vy,
                    const RealT* vz );

  void setResponse( RealT* rx, RealT* ry, RealT* rz );

  const ArrayT<IndexT, 1>& getSortedSurfaceNodeIds() const 
  { 
    return m_sorted_surface_node_ids;
  }

  std::unique_ptr<Viewer> getView() { return std::make_unique<Viewer>(*this); }

private:
  int getDimFromElementType() const;

  template <typename T>
  MultiArrayView<T> createNodalVector( T* x, 
                                        T* y,
                                        T* z ) const;

  ArrayViewT<const IndexT, 2> createConnectivity( IndexT num_elements, 
                                                  const IndexT* connectivity );

  /// sorts unique surface node ids from connectivity and stores them on the mesh object in ascending order
  void sortSurfaceNodeIds();

  IndexT m_mesh_id;                    ///< Mesh Id associated with this data
  InterfaceElementType m_element_type; ///< Type of interface element in mesh
  int m_dim;
  IndexT m_num_nodes;

  MemorySpace m_mem_space; ///< Memory space for mesh data
  int m_allocator_id;      ///< Allocator for mesh data memory
  
  bool m_is_valid;                     ///< True if the mesh is valid

  MeshNodalData m_nodal_fields; ///< method specific nodal fields
  MeshElemData  m_element_data; ///< method/enforcement specific element data

  // Nodal field data
  MultiArrayView<const RealT> m_position; ///< Coordinates of nodes in mesh
  MultiArrayView<const RealT> m_disp;     ///< Nodal displacements
  MultiArrayView<const RealT> m_vel;      ///< Nodal velocity
  MultiArrayView<RealT> m_response;       ///< Nodal responses (forces)

  Array2D<RealT> m_node_n;             ///< Outward unit node normals

  // Element field data
  Array2DView<const IndexT> m_connectivity;  ///< Element connectivity arrays
  Array1D<IndexT> m_sorted_surface_node_ids; ///< List of sorted node ids from connectivity w/o duplicates

  Array2D<RealT> m_c;   ///< Vertex averaged element centroids
  Array2D<RealT> m_n;   ///< Outward unit element normals
  Array1D<RealT> m_face_radius; ///< Face radius used in low level proximity check
  Array1D<RealT> m_area;        ///< Element areas

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

  /*!
  * \brief Computes the face normals and centroids for all faces in the mesh.
  *
  * \param [in] dim Dimension of the problem 
  * \return true if face calculations do not encounter errors or warnings
  * 
  * This routine accounts for warped faces by computing an average normal.
  */
  bool computeFaceData();
  
  /*!
  * \brief Computes average nodal normals for use with mortar methods
  *
  * \note this routine computes average nodal normals for all nodes in the mesh.
  *
  * \param [in] dim Dimension of the problem
  */
  void computeNodalNormals( int const dim );
  
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

  /// Prints information associated with this mesh to \a os
  void print(std::ostream& os) const;
};

//------------------------------------------------------------------------------
template <typename T>
MultiArrayView<T> MeshData::createNodalVector( T* x, T* y, T* z ) const
{
  MultiArrayView<T> host_nodal_vector(m_dim, m_dim);
  host_nodal_vector[0] = ArrayViewT<T>(x, m_num_nodes);
  host_nodal_vector[1] = ArrayViewT<T>(y, m_num_nodes);
  if (m_dim == 3)
  {
    host_nodal_vector[2] = ArrayViewT<T>(z, m_num_nodes);
  }
  return MultiArrayView<T>(host_nodal_vector, m_allocator_id);
}

using MeshManager = DataManager<MeshData>;

} // end namespace tribol

/// \a ostream operator to print a \a MeshData instance to \a os
std::ostream& operator<<(std::ostream& os, const tribol::MeshData& md);

#endif /* SRC_MESH_MESHDATA_HPP_ */
