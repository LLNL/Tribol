// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Tribol Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (MIT)

#include "tribol/mesh/InterfacePairs.hpp"

namespace tribol
{

void InterfacePairs::reserve(integer new_size)
{
   m_pairIndex1.reserve(new_size);
   m_pairIndex2.reserve(new_size);
   m_inContact.reserve(new_size);
}

void InterfacePairs::clear()
{
   m_pairIndex1.clear();
   m_pairIndex2.clear();
   m_inContact.clear();
}


void InterfacePairs::addInterfacePair( InterfacePair const& pair )
{
   // set mesh ids and pair types. This is redundant since these 
   // only need to be set once.
   m_meshId1 = pair.meshId1;
   m_meshId2 = pair.meshId2;
   m_pairType1 = pair.pairType1;
   m_pairType2 = pair.pairType2;

   // add face ids to containers
   m_pairIndex1.push_back(pair.pairIndex1);
   m_pairIndex2.push_back(pair.pairIndex2);

   // set contact boolean container entry
   m_inContact.push_back(pair.inContact);
}

void InterfacePairs::updateInterfacePair( InterfacePair const& pair, 
                                          integer const idx )
{
   m_inContact[ idx ] = pair.inContact;
}

} // namespace tribol


