// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Tribol Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (MIT)

#ifndef SRC_MESH_INTERFACE_PAIRS_HPP_
#define SRC_MESH_INTERFACE_PAIRS_HPP_

#include "tribol/common/ArrayTypes.hpp"
#include "tribol/common/Parameters.hpp"

#include "axom/slic.hpp"

namespace tribol
{

struct InterfacePair
{
   TRIBOL_HOST_DEVICE InterfacePair( int m1, int t1, int i1,
                  int m2, int t2, int i2 )
      : pairId(-1)
      , mesh_id1(m1), pairType1(t1), pairIndex1(i1)
      , mesh_id2(m2), pairType2(t2), pairIndex2(i2) 
      , isContactCandidate(true) {}

   // overload constructor with pair id
   TRIBOL_HOST_DEVICE InterfacePair( int m1, int t1, int i1,
                  int m2, int t2, int i2,
                  int id )
      : pairId(id)
      , mesh_id1(m1), pairType1(t1), pairIndex1(i1)
      , mesh_id2(m2), pairType2(t2), pairIndex2(i2) 
      , isContactCandidate(true) {}

   // overload constructor with boolean indicating contact candidacy
   TRIBOL_HOST_DEVICE InterfacePair( int m1, int t1, int i1,
                  int m2, int t2, int i2,
                  bool isCandidate, int id )
      : pairId(id)
      , mesh_id1(m1), pairType1(t1), pairIndex1(i1)
      , mesh_id2(m2), pairType2(t2), pairIndex2(i2) 
      , isContactCandidate(isCandidate) {}

   // overload constructor to handle zero input arguments
   TRIBOL_HOST_DEVICE InterfacePair() : pairId(-1)
                   , mesh_id1(-1), pairType1(-1), pairIndex1(-1)
                   , mesh_id2(-1), pairType2(-1), pairIndex2(-1)
                   , isContactCandidate(true) {}

   // pair id
   int pairId;

   // mesh id, face type, and face id for face 1
   IndexT mesh_id1;
   int pairType1;
   int pairIndex1;

   // mesh id, face type and face id for face 2
   IndexT mesh_id2;
   int pairType2;
   int pairIndex2;

   // boolean indicating if a binned pair is a contact candidate.
   // A contact candidate is defined as a face-pair that is deemed geometrically proximate 
   // by the binning coarse search, and one that passes the finer computational geometry 
   // (CG) filter/checks. These finer checks identify face-pairs that are intersecting or 
   // nearly intersecting with positive areas of overlap. These checks do not indicate 
   // whether a face-pair is contacting, since the definition of 'contacting' is specific 
   // to a particular contact method and its enforced contstraints. Rather, the CG filter
   // identifies contact candidacy, or face-pairs likely in contact.
   bool isContactCandidate;
};

class InterfacePairs
{
public:
  InterfacePairs()  = default;
  InterfacePairs(int allocator_id);

  void clear();
  void reserve(IndexT capacity);
  void resize(IndexT new_size);

  void addInterfacePair( InterfacePair const& pair );

  int* getPairIndex1Array()
  {
     return m_pairIndex1.data();
  } 

  int* getPairIndex2Array()
  {
     return m_pairIndex2.data();
  } 

  bool* getContactArray()
  {
     return m_isContactCandidate.data();
  } 

  TRIBOL_HOST_DEVICE void updateInterfacePair( InterfacePair const& pair,
                            int const idx );

  TRIBOL_HOST_DEVICE void setMeshId( int side, int id )
  {
     SLIC_ERROR_IF( side !=1 && side !=2, "mesh side not 1 or 2.");
     if (side == 1) m_mesh_id1 = id;
     if (side == 2) m_mesh_id2 = id;
  }

  TRIBOL_HOST_DEVICE void setPairType( int side, int type )
  {
     SLIC_ERROR_IF( side !=1 && side !=2, "mesh side not 1 or 2.");
     if (side == 1) m_pairType1 = type;
     if (side == 2) m_pairType2 = type;
  }

  TRIBOL_HOST_DEVICE InterfacePair getInterfacePair(IndexT idx) const;

  TRIBOL_HOST_DEVICE IndexT getNumPairs() const { return m_pairIndex1.size(); }

private:

  // A list of interface pairs is defined to be between 
  // two meshes with the same elements type on each mesh
  IndexT m_mesh_id1;
  IndexT m_mesh_id2;
  int m_pairType1;
  int m_pairType2;

  int m_allocator_id;

  ArrayT<int> m_pairIndex1;
  ArrayT<int> m_pairIndex2;
  ArrayT<bool> m_isContactCandidate;
};

} /* namespace tribol */

#endif /* SRC_MESH_INTERFACE_PAIRS_HPP_ */
