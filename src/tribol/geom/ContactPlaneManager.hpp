// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Tribol Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (MIT)

#ifndef SRC_GEOM_CONTACTPLANEMANAGER_HPP_
#define SRC_GEOM_CONTACTPLANEMANAGER_HPP_

#include "tribol/types.hpp"
#include "tribol/geom/ContactPlane.hpp"

namespace tribol
{

class ContactPlaneManager
{
protected:

  int m_spaceDim; ///< Space dimension

  IndexT m_numContactPlanes; ///< Total number of contact planes
  IndexT m_allocatedNumContactPlanes; ///< Total number of allocated contact planes

  bool m_resized;

public:

  ArrayT<IndexT> m_numFaces; ///< Number of constituent faces

  ArrayT<IndexT> m_meshId1; ///< Mesh id for first face
  ArrayT<IndexT> m_meshId2; ///< Mesh id for second face

  ArrayT<IndexT> m_fId1; ///< First face id
  ArrayT<IndexT> m_fId2; ///< Second face id

  ArrayT<RealT> m_pressure; ///< Contact pressure
  ArrayT<RealT> m_ratePressure; ///< Velocity gap pressure

  ArrayT<bool> m_inContact; ///< True if faces are in contact

  ArrayT<bool> m_interpenOverlap; ///< True if interpenetration overlap calculation

  ArrayT<RealT> m_cX; ///< Global x-coordinates of contact plane point
  ArrayT<RealT> m_cY; ///< Global y-coordinates of contact plane point
  ArrayT<RealT> m_cZ; ///< Global z-coordinates of contact plane point

  ArrayT<RealT> m_e1X; ///< Global x-component of first in-plane basis vector
  ArrayT<RealT> m_e1Y; ///< Global y-component of first in-plane basis vector
  ArrayT<RealT> m_e1Z; ///< Global z-component of first in-plane basis vector

  ArrayT<RealT> m_e2X; ///< Global x-component of second in-plane basis vector
  ArrayT<RealT> m_e2Y; ///< Global y-component of second in-plane basis vector
  ArrayT<RealT> m_e2Z; ///< Global z-component of second in-plane basis vector

  ArrayT<RealT> m_nX; ///< Global x-component of contact plane unit normal
  ArrayT<RealT> m_nY; ///< Global y-component of contact plane unit normal
  ArrayT<RealT> m_nZ; ///< Global z-component of contact plane unit normal

  ArrayT<IndexT> m_numPolyVert; ///< Number of overlap polygon vertices

  ArrayT<RealT*> m_polyLocX; ///< Local x-coordinates of overlap polygon's vertices
  ArrayT<RealT*> m_polyLocY; ///< Local y-coordinates of overlap polygon's vertices

  ArrayT<RealT*> m_polyX; ///< Global x-coordinates of overlap polygon's vertices
  ArrayT<RealT*> m_polyY; ///< Global y-coordinates of overlap polygon's vertices
  ArrayT<RealT*> m_polyZ; ///< Global z-coordinates of overlap polygon's vertices

  ArrayT<RealT> m_overlapCX; ///< Local x-coordinates of overlap centroid
  ArrayT<RealT> m_overlapCY; ///< Local y-coordinates of overlap centroid

  ArrayT<RealT> m_gap; ///< Contact plane kinematic gap
  ArrayT<RealT> m_velGap; ///< Velocity (rate) gap
  ArrayT<RealT> m_gapTol; ///< Gap tolerance

  ArrayT<IndexT> m_numInterpenPoly1Vert; ///< Number of vertices in face 1 interpen. polygon
  ArrayT<RealT*> m_interpenPoly1X; ///< Local x-coordinates of face 1 interpenetrating polygon 
  ArrayT<RealT*> m_interpenPoly1Y; ///< Local y-coordinates of face 1 interpenetrating polygon 

  ArrayT<IndexT> m_numInterpenPoly2Vert; ///< Number of vertices in face 2 interpen. polygon
  ArrayT<RealT*> m_interpenPoly2X; ///< Local x-coordinates of face 2 interpenetrating polygon
  ArrayT<RealT*> m_interpenPoly2Y; ///< Local y-coordinates of face 2 interpenetrating polygon

  ArrayT<RealT*> m_interpenG1X; ///< Global x-coordinate of face 1 interpenetrating polygon
  ArrayT<RealT*> m_interpenG1Y; ///< Global y-coordinate of face 1 interpenetrating polygon
  ArrayT<RealT*> m_interpenG1Z; ///< Global z-coordinate of face 1 interpenetrating polygon

  ArrayT<RealT*> m_interpenG2X; ///< Global x-coordinate of face 2 interpenetrating polygon
  ArrayT<RealT*> m_interpenG2Y; ///< Global y-coordinate of face 2 interpenetrating polygon
  ArrayT<RealT*> m_interpenG2Z; ///< Global z-coordiante of face 2 interpenetrating polygon

  ArrayT<RealT> m_areaFrac; ///< Area fraction for polygon overlap inclusion
  ArrayT<RealT> m_areaMin; ///< Minimum overlap area to include in the active set
  ArrayT<RealT> m_area; ///< Overlap area
  ArrayT<RealT> m_interpenArea; ///< Area of the interpenetrating overlap

  ArrayT<RealT> m_cXf1; ///< Global x-coord of the contact plane centroid projected to face 1
  ArrayT<RealT> m_cYf1; ///< Global y-coord of the contact plane centroid projected to face 1
  ArrayT<RealT> m_cZf1; ///< Global z-coord of the contact plane centroid projected to face 1

  ArrayT<RealT> m_cXf2; ///< Global x-coord of the contact plane centroid projected to face 2
  ArrayT<RealT> m_cYf2; ///< Global y-coord of the contact plane centroid projected to face 2
  ArrayT<RealT> m_cZf2; ///< Global z-coord of the contact plane centroid projected to face 2

  ArrayT<RealT*> m_segX; ///< Global x-coordinate of the overlap vertices (2D)
  ArrayT<RealT*> m_segY; ///< Global y_coordinate of the overlap vertices (2D)

public:
  /*!
  * \brief Destructor
  *
  */
  ~ContactPlaneManager() { if( m_numContactPlanes > 0 ) clearCPManager(); } ;

  /*!
  * \brief Returns a reference to the ContactPlaneManager instance.
  * \return C a reference to the ContactPlaneManager instance.
  *
  */
  static ContactPlaneManager & getInstance(); 

  /*!
  * /brief Returns the number of contact planes
  *
  * \return Number of contact planes 
  */
  int size() const { return m_numContactPlanes; }

  /*!
   * \brief Wrapper routine for adding 2D and 3D contact planes to the list.
   *
   * \param [in] newSize New size of the container
   */
  void resize( IndexT const newSize );

  /*!
   * \brief Wrapper routine to allocate space for newSize 2D or 
   * 3D contact planes.
   *
   * \param [in] newSize New size of the container
   */
  void reserve( IndexT const newSize );

  /*!
  * \brief 3D routine to allocate space for newSize contact planes
  *
  * \param [in] newSize New size of the container
  */
  void reserve3D( IndexT const newSize );

  /*!
  * \brief 2D routine to allocate space for newSize contact planes
  *
  * \param [in] newSize New size of the container
  */
  void reserve2D( IndexT const newSize );

  /*!
  * \brief Set the space dimension
  *
  * \param [in] dimension Space dimension
  */
  void setSpaceDim( int dimension) { m_spaceDim = dimension; };

  /*!
  * \brief Get the space dimension
  *
  * \return Space dimension
  */
  int getSpaceDim() const { return m_spaceDim; };

  /*!
  * \brief Copies 3D contact plane data over to manager
  *
  * \param [in] cp Contact plane object
  */
  void setContactPlaneData( const ContactPlane3D& cp );

  /*!
  * \brief Adds 3D contact plane data to manager
  *
  * \param [in] cp 3D contact plane object
  */
  void addContactPlane( const ContactPlane3D& cp );

  /*!
  * \brief Copies 2D contact plane data over to manager
  *
  * \param [in] cp Contact plane object
  */
  void setContactPlaneData( const ContactPlane2D& cp );

  /*!
  * \brief Adds 2D contact plane data to manager
  *
  * \param [in] cp 2D contact plane object
  */
  void addContactPlane( const ContactPlane2D& cp );

  /*!
  * \brief populates pointer to array of stacked (x,y,z) coordinates of 
  *        contact plane overlap vertices 
  *
  * \param [in] id contact plane id
  * \param [in] numVerts number of overlap vertices
  * \param [in/out] coords pointer to array of stacked (x,y,z) coordinates of overlap vertices 
  *
  * \pre coords points to a 1D array of size (problem dimension) x numVerts.
  */
  void getContactPlaneOverlapVerts( int const id, int const numVerts, RealT * coords );

  /*!
  * \brief returns pointer to array of stacked components of the contact 
  *        plane normal vector
  *
  * \param [in] id contact plane id
  * \param [in] dim dimension of the problem
  * \param [in/out] coords pointer to array of stacked normal components 
  *
  */
  void getContactPlaneNormal( int const id, int const dim, RealT * nrml );

  /*!
  * \brief populates pointer to array of stacked (x,y,z) coordinates of 
  *        face vertices as projected onto the contact plane
  *
  * \param [in] id contact plane id
  * \param [in] faceId local face integer id
  * \param [in/out] coords pointer to array of stacked (x,y,z) coordinates of projected face vertices 
  *
  * \pre coords points to a 1D array of size (problem dimension) x number of face vertices
  * \pre faceId is 0 or 1 for first or second face, respectively
  */
  void getProjectedFaceCoords( int const id, int const faceId,
                               RealT * coords );

  /*!
  * \brief clear contact plane manager 
  *
  * \note this clears data and deallocates/nullptr for allocatable arrays
  *
  */
  void clearCPManager();

private:
  /*!
  * \brief 2D resize routine
  *
  * \note Assumes the new size \a m_numContactPlanes has already been set
  */
  void resize2D();

  /*!
  * \brief 3D resize routine
  *
  * \note Assumes the new size \a m_numContactPlanes has already been set
  */
  void resize3D();

private:

  /*!
  * \brief Private constructor 
  *
  */
  ContactPlaneManager();

};

}
#endif /* SRC GEOM_CONTACTPLANEMANAGER_HPP_ */
