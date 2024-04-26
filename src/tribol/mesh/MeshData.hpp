// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Tribol Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (MIT)

#ifndef SRC_MESH_MESHDATA_HPP_
#define SRC_MESH_MESHDATA_HPP_

#include "tribol/common/ArrayTypes.hpp"
#include "tribol/common/LoopExec.hpp"
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
  int m_num_nodes;

  /////////////////////////
  // MORTAR NODAL FIELDS //
  /////////////////////////
  ArrayViewT<RealT>       m_node_gap;      ///< scalar nodal gap (used on nonmortar mesh) 
  ArrayViewT<const RealT> m_node_pressure; ///< scalar nodal pressure (used on nonmortar mesh) 

  bool m_is_gap_computed      {false}; ///< true if the nodal gaps have been computed
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
  virtual MemorySpace getMemorySpace() const = 0;

  IndexT meshId() { return m_mesh_id; }

  InterfaceElementType getElementType() const { return m_element_type; }

  int dimension() const { return m_dim; }

  bool& isMeshValid() { return m_is_valid; }

  IndexT numberOfNodes() const { return m_num_nodes; };

  IndexT numberOfElements() const { return m_num_elements; };

  TRIBOL_HOST_DEVICE virtual IndexT numberOfNodesPerElement() const = 0;

  TRIBOL_HOST_DEVICE virtual IndexT getGlobalNodeId(IndexT element_id, IndexT local_node_id) const = 0;

  virtual bool hasVelocity() const = 0;

  virtual bool hasNodalNormals() const = 0;

  virtual const RealT* getPosition( int dim ) const = 0;

  virtual const RealT* getDisplacement( int dim ) const = 0;

  virtual const RealT* getVelocity( int dim ) const = 0;

  virtual RealT* getResponse( int dim ) const = 0;

  virtual void setPosition( const RealT* x,
                            const RealT* y,
                            const RealT* z ) = 0;

  virtual void setDisplacement( const RealT* ux,
                                const RealT* uy,
                                const RealT* uz ) = 0;

  virtual void setVelocity( const RealT* ux,
                            const RealT* uy,
                            const RealT* uz ) = 0;
                            
  virtual void setResponse( RealT* rx, RealT* ry, RealT* rz ) = 0;

  virtual const IndexT* getConnectivityData() const = 0;
  
  virtual const RealT* getNodalNormals(IndexT dim) const = 0;

  virtual const IndexT* getSortedSurfaceNodeIdsData() const = 0;

  virtual IndexT getSortedSurfaceNodeIdsSize() const = 0;

  virtual const RealT* getElementNormals(IndexT dim) const = 0;

  virtual const RealT* getElementCentroids(IndexT dim) const = 0;

  virtual const RealT* getFaceRadiusData() const = 0;

  virtual const RealT* getElementAreaData() const = 0;

  MeshNodalData& getNodalFields() { return m_nodal_fields; }

  MeshElemData& getElementData() { return m_element_data; }

  const MeshElemData& getElementData() const { return m_element_data; }

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
  virtual bool computeFaceData() = 0;
  
  /*!
  * \brief Computes average nodal normals for use with mortar methods
  *
  * \note this routine computes average nodal normals for all nodes in the mesh.
  *
  * \param [in] dim Dimension of the problem
  */
  virtual void computeNodalNormals( int const dim ) = 0;
  
  /*!
  *
  * \brief returns pointer to array of stacked nodal coordinates for given face
  *
  * \param [in/out] coords pointer to an array of stacked (x,y,z) nodal coordinates
  *
  */
  virtual void getFaceCoords( int const faceId, RealT * coords ) = 0;

  /*!
  *
  * \brief returns pointer to array of stacked nodal velocities for given face
  *
  * \param [in/out] nodalVel pointer to an array of stacked (x,y,z) nodal velocities
  *
  */
  virtual void getFaceNodalVelocities( int const faceId, RealT * nodalVel ) = 0;

  /*!
  *
  * \brief returns pointer to array of stacked normal components 
  *
  * \param [in] faceId integer id of face
  * \param [in] dim dimension of the problem
  * \param [in/out] nrml pointer to array of stacked components of the face normal vector
  *
  */
  virtual void getFaceNormal( int const faceId, int const dim, RealT * nrml ) = 0;

  virtual void print(std::ostream& os) const = 0;

protected:
  MeshData( IndexT mesh_id, InterfaceElementType element_type,
    IndexT num_nodes, IndexT num_elements );

  int getDimFromElementType();

  IndexT m_mesh_id;                    ///< Mesh Id associated with this data
  InterfaceElementType m_element_type; ///< Type of interface element in mesh
  int m_dim;                           ///< Mesh dimension

  IndexT m_num_nodes;                  ///< Number of nodes in the mesh
  IndexT m_num_elements;               ///< Number of elements in the mesh
  
  bool m_is_valid;                     ///< True if the mesh is valid

  MeshNodalData m_nodal_fields; ///< method specific nodal fields
  MeshElemData  m_element_data; ///< method/enforcement specific element data
};

template <MemorySpace MSPACE = MemorySpace::Host>
class MeshDataBySpace : public MeshData
{
public:
  /*!
  * \brief Constructor 
  *
  */
  MeshDataBySpace( IndexT mesh_id, IndexT num_elements, IndexT num_nodes,
                   const IndexT* connectivity, InterfaceElementType element_type,
                   const RealT* x, const RealT* y, const RealT* z );

  MemorySpace getMemorySpace() const override { return MSPACE; }

  bool hasVelocity() const override { return !m_vel.empty(); }

  bool hasNodalNormals() const override { return !m_node_n.empty(); }

  const RealT* getPosition(IndexT dim) const override
  {
    return dim >= m_position.size() ? nullptr : m_position[dim].data();
  }

  const VectorArrayView<const RealT, MSPACE>& getPosition() const 
  {
    return m_position;
  }

  const RealT* getDisplacement(IndexT dim) const override
  {
    return dim >= m_disp.size() ? nullptr : m_disp[dim].data();
  }

  const VectorArrayView<const RealT, MSPACE>& getDisplacement() const
  { 
    return m_disp;
  }

  const RealT* getVelocity(IndexT dim) const override
  {
    return dim >= m_vel.size() ? nullptr : m_vel[dim].data();
  }

  const VectorArrayView<const RealT, MSPACE>& getVelocity() const
  { 
    return m_vel;
  }

  RealT* getResponse(IndexT dim) const override
  {
    return dim >= m_response.size() ? nullptr : m_response[dim].data();
  }

  const VectorArrayView<RealT, MSPACE>& getResponse() const
  {
    return m_response;
  }

  VectorArrayView<RealT, MSPACE>& getResponse()
  {
    return m_response;
  }

  const RealT* getNodalNormals(IndexT dim) const override
  {
    return dim >= m_node_n.size() ? nullptr : m_node_n[dim].data();
  }

  const VectorArray<RealT, MSPACE>& getNodalNormals() const { return m_node_n; }

  void setPosition( const RealT* x,
                    const RealT* y,
                    const RealT* z ) override;

  void setDisplacement( const RealT* ux,
                        const RealT* uy,
                        const RealT* uz ) override;

  void setVelocity( const RealT* vx,
                    const RealT* vy,
                    const RealT* vz ) override;

  void setResponse( RealT* rx, RealT* ry, RealT* rz ) override;

  const IndexT* getConnectivityData() const override { return m_connectivity.data(); }

  const ArrayViewT<const IndexT, 2, MSPACE>& getConnectivity() const
  {
    return m_connectivity;
  }

  const IndexT* getSortedSurfaceNodeIdsData() const override
  {
    return m_sorted_surface_node_ids.data();
  }

  IndexT getSortedSurfaceNodeIdsSize() const override
  {
    return m_sorted_surface_node_ids.size();
  }

  const ArrayT<IndexT, 1, MSPACE>& getSortedSurfaceNodeIds() const 
  { 
    return m_sorted_surface_node_ids;
  }

  TRIBOL_HOST_DEVICE IndexT numberOfNodesPerElement() const override
  { 
    return m_connectivity.shape()[1];
  }

  TRIBOL_HOST_DEVICE IndexT getGlobalNodeId(IndexT element_id, IndexT local_node_id) const override
  {
    return m_connectivity(element_id, local_node_id);
  }

  const RealT* getElementNormals(IndexT dim) const override
  {
    return dim >= m_n.size() ? nullptr : m_n[dim].data();
  }

  const VectorArray<RealT, MSPACE>& getElementNormals() const { return m_n; }

  const RealT* getElementCentroids(IndexT dim) const override
  {
    return dim >= m_c.size() ? nullptr : m_c[dim].data();
  }

  const VectorArray<RealT, MSPACE>& getElementCentroids() const { return m_c; }

  const RealT* getFaceRadiusData() const override { return m_face_radius.data(); }

  const ScalarArray<RealT, MSPACE>& getFaceRadius() const { return m_face_radius; }

  const RealT* getElementAreaData() const override { return m_area.data(); }

  const ScalarArray<RealT, MSPACE>& getElementAreas() const { return m_area; }

private:
  ArrayViewT<const IndexT, 2, MSPACE> createConnectivity( IndexT num_elements, 
                                                          const IndexT* connectivity );

  /// sorts unique surface node ids from connectivity and stores them on the mesh object in ascending order
  void sortSurfaceNodeIds();

  // Nodal field data
  VectorArrayView<const RealT, MSPACE> m_position; ///< Coordinates of nodes in mesh
  VectorArrayView<const RealT, MSPACE> m_disp;     ///< Nodal displacements
  VectorArrayView<const RealT, MSPACE> m_vel;      ///< Nodal velocity

  VectorArrayView<RealT, MSPACE> m_response;       ///< Nodal responses (forces)
  VectorArray<RealT, MSPACE> m_node_n;             ///< Outward unit node normals

  // Element field data
  ArrayViewT<const IndexT, 2, MSPACE> m_connectivity;  ///< Element connectivity arrays
  ArrayT<IndexT, 1, MSPACE> m_sorted_surface_node_ids; ///< List of sorted node ids from connectivity w/o duplicates

  VectorArray<RealT, MSPACE> m_n;   ///< Outward unit element normals
  VectorArray<RealT, MSPACE> m_c;   ///< Vertex averaged element centroids
  ScalarArray<RealT, MSPACE> m_face_radius; ///< Face radius used in low level proximity check
  ScalarArray<RealT, MSPACE> m_area;        ///< Element areas

public:
  template <typename T>
  VectorArrayView<T, MSPACE> createNodalVector( T* x, 
                                                T* y,
                                                T* z ) const;

  /*!
  * \brief Computes the face normals and centroids for all faces in the mesh.
  *
  * \param [in] dim Dimension of the problem 
  * \return true if face calculations do not encounter errors or warnings
  * 
  * This routine accounts for warped faces by computing an average normal.
  */
  bool computeFaceData() override;
  
  /*!
  * \brief Computes average nodal normals for use with mortar methods
  *
  * \note this routine computes average nodal normals for all nodes in the mesh.
  *
  * \param [in] dim Dimension of the problem
  */
  void computeNodalNormals( int const dim ) override;
  
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
  TRIBOL_HOST_DEVICE RealT computeFaceRadius( int faceId );
  
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
  void getFaceCoords( int const faceId, RealT * coords ) override;

  /*!
  *
  * \brief returns pointer to array of stacked nodal velocities for given face
  *
  * \param [in/out] nodalVel pointer to an array of stacked (x,y,z) nodal velocities
  *
  */
  void getFaceNodalVelocities( int const faceId, RealT * nodalVel ) override;

  /*!
  *
  * \brief returns pointer to array of stacked normal components 
  *
  * \param [in] faceId integer id of face
  * \param [in] dim dimension of the problem
  * \param [in/out] nrml pointer to array of stacked components of the face normal vector
  *
  */
  void getFaceNormal( int const faceId, int const dim, RealT * nrml ) override;

  /// Prints information associated with this mesh to \a os
  void print(std::ostream& os) const override;
};

//------------------------------------------------------------------------------
template <MemorySpace MSPACE>
template <typename T>
VectorArrayView<T, MSPACE> MeshDataBySpace<MSPACE>::createNodalVector( T* x,
                                                                       T* y,
                                                                       T* z ) const
{
  VectorArrayView<T, MSPACE> nodal_vector(m_dim, m_dim);
  // auto nodal_vector_data = nodal_vector.data();
  // T* vector_temp[3]{x, y, z};

  // forAllExec<toExecutionMode<MSPACE>::value>(
  //   m_dim,
  //   [this, nodal_vector_data, vector_temp] TRIBOL_HOST_DEVICE (const IndexT i) {
  //     axom::StackArray<IndexT, 1> num_nodes_array;
  //     num_nodes_array[0] = numberOfNodes();
  //     nodal_vector_data[i] = ArrayViewT<T, 1, MSPACE>(vector_temp[i], num_nodes_array);
  //   }
  // );
  ArrayT<ArrayViewT<T, 1, MSPACE>, 1, MemorySpace::Dynamic> host_nodal_vector(m_dim, m_dim);
  host_nodal_vector[0] = ArrayViewT<T, 1, MSPACE>(x, numberOfNodes());
  host_nodal_vector[1] = ArrayViewT<T, 1, MSPACE>(y, numberOfNodes());
  if (m_dim == 3)
  {
    host_nodal_vector[2] = ArrayViewT<T, 1, MSPACE>(z, numberOfNodes());
  }
  nodal_vector = std::move(host_nodal_vector);

  return nodal_vector;
}

// use a unique_ptr for polymorphism
using MeshManager = DataManager<std::unique_ptr<MeshData>>;

} // end namespace tribol

/// \a ostream operator to print a \a MeshData instance to \a os
std::ostream& operator<<(std::ostream& os, const tribol::MeshData& md);

#endif /* SRC_MESH_MESHDATA_HPP_ */
