// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Tribol Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (MIT)

#include "tribol/mesh/InterfacePairs.hpp"
#include "tribol/common/ExecModel.hpp"

namespace tribol
{

InterfacePairs::InterfacePairs(MemorySpace mem_space)
: m_allocator_id(getResourceAllocatorID(mem_space)),
  m_pairIndex1(0, 0, m_allocator_id),
  m_pairIndex2(0, 0, m_allocator_id),
  m_isContactCandidate(0, 0, m_allocator_id)
{}

void InterfacePairs::reserve(int new_size)
{
   m_pairIndex1.reserve(new_size);
   m_pairIndex2.reserve(new_size);
   m_isContactCandidate.reserve(new_size);
}

void InterfacePairs::clear()
{
   m_pairIndex1.clear();
   m_pairIndex2.clear();
   m_isContactCandidate.clear();
}


TRIBOL_HOST_DEVICE void InterfacePairs::addInterfacePair( InterfacePair const& pair )
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

TRIBOL_HOST_DEVICE void InterfacePairs::updateInterfacePair( InterfacePair const& pair, 
                                          int const idx )
{
   m_isContactCandidate[ idx ] = pair.isContactCandidate;
}

TRIBOL_HOST_DEVICE InterfacePair InterfacePairs::getInterfacePair(IndexT idx) const
{
   return InterfacePair {
      m_mesh_id1, m_pairType1, m_pairIndex1[idx],
      m_mesh_id2, m_pairType2, m_pairIndex2[idx], 
      m_isContactCandidate[idx], idx };
}

} // namespace tribol


