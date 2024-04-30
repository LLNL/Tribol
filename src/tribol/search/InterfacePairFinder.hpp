/*
 *******************************************************************************
 * Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 *******************************************************************************
 */
#ifndef TRIBOL_SEARCH_INTERFACE_PAIR_FINDER_HPP_
#define TRIBOL_SEARCH_INTERFACE_PAIR_FINDER_HPP_

#include "tribol/common/Parameters.hpp"
#include "tribol/mesh/MeshData.hpp"

namespace tribol
{

// Forward Declarations
class CouplingScheme;
struct InterfacePair;
class SearchBase;


/// Free functions

/*!
 * \brief Basic geometry/proximity checks for face pairs
 *
 * \param [in] pairIndex1 index of 1st element in pair
 * \param [in] pairIndex2 index of 2nd element in pair
 * \param [in] meshId1 mesh ID for 1st element in pair
 * \param [in] meshId2 mesh ID for 2nd element in pair
 * \param [in] pMesh1 pointer to data for mesh with ID meshId1
 * \param [in] pMesh2 pointer to data for mesh with ID meshId2
 * \param [in] mode ContactMode
 *
 */
//bool geomFilter( InterfacePair const & iPair, ContactMode const mode );
TRIBOL_HOST_DEVICE bool geomFilter( const IndexT pairIndex1, const IndexT pairIndex2,
                                    const IndexT meshId1, const IndexT meshId2, 
                                    const MeshData* const pMesh1, const MeshData* const pMesh2,
                                    ContactMode const mode );

/*!
 * \class InterfacePairFinder
 *
 * \brief This class finds pairs of interfering elements in the meshes
 * referred to by the CouplingScheme
 */
class InterfacePairFinder
{
public:
   InterfacePairFinder(CouplingScheme* cs);

   ~InterfacePairFinder();

   /*!
    * Initializes structures for the candidate search
    */
   void initialize();

   /*!
    * Computes the interacting interface pairs between the meshes
    * specified in \a m_couplingScheme
    */
   void findInterfacePairs();

private:

   CouplingScheme* m_couplingScheme;
   SearchBase* m_search;  // The search strategy
};

} // end namespace tribol



#endif /* TRIBOL_SEARCH_INTERFACE_PAIR_FINDER_HPP_ */
