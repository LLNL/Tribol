// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Tribol Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (MIT)

#include "tribol/mesh/InterfacePairs.hpp"
#include "tribol/common/ExecModel.hpp"

namespace tribol
{

InterfacePairs::InterfacePairs(int allocator_id)
  : m_allocator_id( allocator_id )
  , m_pairIndex1( 0, 0, m_allocator_id )
  , m_pairIndex2( 0, 0, m_allocator_id )
  , m_isContactCandidate( 0, 0, m_allocator_id )
{}

void InterfacePairs::reserve(IndexT capacity)
{
  m_pairIndex1.reserve(capacity);
  m_pairIndex2.reserve(capacity);
  m_isContactCandidate.reserve(capacity);
}

void InterfacePairs::resize(IndexT new_size)
{
   m_pairIndex1.resize(new_size);
   m_pairIndex2.resize(new_size);
   m_isContactCandidate.resize(new_size);
}

void InterfacePairs::clear()
{
   m_pairIndex1.clear();
   m_pairIndex2.clear();
   m_isContactCandidate.clear();
}

void InterfacePairs::addInterfacePair( InterfacePair const& pair )
{
   // set mesh ids and pair types. This is redundant since these 
   // only need to be set once.
  //  m_mesh_id1 = pair.mesh_id1;
  //  m_mesh_id2 = pair.mesh_id2;
  //  m_pairType1 = pair.pairType1;
  //  m_pairType2 = pair.pairType2;

   // add face ids to containers
   m_pairIndex1.resize(m_pairIndex1.size() + 1, pair.pairIndex1);
   m_pairIndex1.push_back(pair.pairIndex1);
   m_pairIndex2.push_back(pair.pairIndex2);

   // set contact candidate boolean container entry
   m_isContactCandidate.push_back(pair.isContactCandidate);
}

InterfacePairs::Viewer::Viewer(const InterfacePairs& pairs)
  : m_mesh_id1( pairs.m_mesh_id1 )
  , m_mesh_id2( pairs.m_mesh_id2 )
  , m_pairType1( pairs.m_pairType1 )
  , m_pairType2( pairs.m_pairType2 )
  , m_pairIndex1( pairs.m_pairIndex1 )
  , m_pairIndex2( pairs.m_pairIndex2 )
  , m_isContactCandidate( pairs.m_isContactCandidate )
{}

TRIBOL_HOST_DEVICE InterfacePair InterfacePairs::Viewer::getInterfacePair(IndexT idx) const
{
   return {
      m_mesh_id1, m_pairType1, m_pairIndex1[idx],
      m_mesh_id2, m_pairType2, m_pairIndex2[idx], 
      m_isContactCandidate[idx], idx };
}

TRIBOL_HOST_DEVICE void InterfacePairs::Viewer::updateInterfacePair(
  InterfacePair const& pair, int const idx ) const
{
   m_isContactCandidate[ idx ] = pair.isContactCandidate;
}

} // namespace tribol


