// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Tribol Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (MIT)

#include "tribol/search/InterfacePairFinder.hpp"

#include "tribol/types.hpp"
#include "tribol/mesh/CouplingScheme.hpp"
#include "tribol/mesh/MeshData.hpp"
#include "tribol/mesh/InterfacePairs.hpp"
#include "tribol/common/Parameters.hpp"
#include "tribol/utils/Math.hpp"

#include "axom/slic.hpp"
#include "axom/slam.hpp"
#include "axom/primal.hpp"
#include "axom/spin.hpp"

// Define some namespace aliases to help with axom usage
namespace slam = axom::slam;
namespace primal = axom::primal;
namespace spin = axom::spin;

namespace tribol
{

/*!
 *  Perform geometry/proximity checks 1-4
 */
bool geomFilter( InterfacePair & iPair, ContactMode const mode )
{
   // alias variables off the InterfacePair
   const IndexT& mesh_id1 = iPair.mesh_id1;
   const IndexT& mesh_id2 = iPair.mesh_id2;
   int const & faceId1 = iPair.pairIndex1;
   int const & faceId2 = iPair.pairIndex2;

   /// CHECK #1: Check to make sure the two face ids are not the same
   ///           and the two mesh ids are not the same.
   if ((mesh_id1 == mesh_id2) && (faceId1 == faceId2))
   {
      iPair.isContactCandidate = false;
      return iPair.isContactCandidate;
   }

   // get instance of mesh manager
   MeshManager & meshManager = MeshManager::getInstance();

   // get instance of mesh data
   MeshData& mesh1 = *meshManager.at(mesh_id1);
   MeshData& mesh2 = *meshManager.at(mesh_id2);
   int dim = mesh1.dimension();

   /// CHECK #2: Check to make sure faces don't share a common
   ///           node for the case where mesh_id1 = mesh_id2.
   ///           We want to preclude two adjacent faces from interacting.
   if (mesh_id1 == mesh_id2)
   {
      for (int i=0; i<mesh1.numberOfNodesPerElement(); ++i)
      {
         int node1 = mesh1.getGlobalNodeId(faceId1, i);
         for (int j=0; j<mesh2.numberOfNodesPerElement(); ++j)
         {
            int node2 = mesh2.getGlobalNodeId(faceId2, j);
            if (node1 == node2)
            {
              iPair.isContactCandidate = false;
              return iPair.isContactCandidate;
            }
         }
      }
   }

   /// CHECK #3: Check that face normals are opposing up to some tolerance.
   ///           This uses a hard coded normal tolerance for this check.
   RealT nrmlTol = -0.173648177; // taken as cos(100) between face pair

   RealT m_nZ1, m_nZ2;
   if ( dim == 3 )
   {
      m_nZ1 = mesh1.getElementNormals(2)[ faceId1 ];
      m_nZ2 = mesh2.getElementNormals(2)[ faceId2 ];
   }
   else
   {
      m_nZ1 = 0.;
      m_nZ2 = 0.;
   }

   RealT nrmlCheck = mesh1.getElementNormals(0)[faceId1] * mesh2.getElementNormals(0)[faceId2] +
                    mesh1.getElementNormals(1)[faceId1] * mesh2.getElementNormals(1)[faceId2] +
                    m_nZ1 * m_nZ2;

   // check normal projection against tolerance
   if (nrmlCheck > nrmlTol) {
      SLIC_DEBUG("Face does not pass normal tolerance");
      iPair.isContactCandidate = false;
      return iPair.isContactCandidate;
   }

   // TODO this may still add faces on opposing sides of thin-walled structures when using
   // auto-contact, which is not correct. Consider comparing against element thicknesses; thus,
   // requiring element thickness for auto contact
   //
   /// CHECK #4 (3D): Perform radius check, which involves seeing if
   ///                the distance between the two face vertex averaged
   ///                centroid is less than the sum of the two face radii.
   ///                The face radii are taken to be the magnitude of the
   ///                longest vector from that face's vertex averaged
   ///                centroid to one its nodes.
   RealT offset_tol = 0.05;
   if (dim == 3)
   {
      RealT r1 = mesh1.getFaceRadiusData()[ faceId1 ];
      RealT r2 = mesh2.getFaceRadiusData()[ faceId2 ];

      // set maximum offset of face centroids for inclusion
      RealT distMax = r1 + r2; // default is sum of face radii

      // check if the contact mode is conforming, in which case the
      // faces are supposed to be aligned
      if (mode == SURFACE_TO_SURFACE_CONFORMING)
      {
         // use 5% of max face radius for conforming case as
         // tolerance on face offsets
         distMax *= offset_tol;
      }

      // compute the distance between the two face centroids
      RealT distX = mesh2.getElementCentroids(0)[ faceId2 ] - mesh1.getElementCentroids(0)[ faceId1 ];
      RealT distY = mesh2.getElementCentroids(1)[ faceId2 ] - mesh1.getElementCentroids(1)[ faceId1 ];
      RealT distZ = mesh2.getElementCentroids(2)[ faceId2 ] - mesh1.getElementCentroids(2)[ faceId1 ];

      RealT distMag = magnitude(distX, distY, distZ );

      if (distMag >= (distMax)) {
         iPair.isContactCandidate = false;
         return iPair.isContactCandidate;
      }
   } // end of dim == 3
   else if (dim == 2)
   {
      // get 1/2 edge length off the mesh data
      RealT e1 = 0.5 * mesh1.getElementAreaData()[ faceId1 ];
      RealT e2 = 0.5 * mesh2.getElementAreaData()[ faceId2 ];

      // set maximum offset of edge centroids for inclusion
      RealT distMax = e1 + e2; // default is sum of 1/2 edge lengths

      // check if the contact mode is conforming, in which case the
      // edges are supposed to be aligned
      if (mode == SURFACE_TO_SURFACE_CONFORMING)
      {
         // use 5% of max face radius for conforming case as
         // tolerance on face offsets
         distMax *= offset_tol;
      }

      // compute the distance between the two edge centroids
      RealT distX = mesh2.getElementCentroids(0)[ faceId2 ] - mesh1.getElementCentroids(0)[ faceId1 ];
      RealT distY = mesh2.getElementCentroids(1)[ faceId2 ] - mesh1.getElementCentroids(1)[ faceId1 ];

      RealT distMag = magnitude(distX, distY);

      if (distMag >= (distMax))
      {
         iPair.isContactCandidate = false;
         return iPair.isContactCandidate;
      }
   } // end of dim == 2

   // if we made it here we passed all checks
   iPair.isContactCandidate = true;
   return iPair.isContactCandidate;

} // end geomFilter()


/*!
 * Wraps a MeshData instance to simplify operations like accessing
 * vertex positions and element bounding boxes
 *
 * \tparam D The spatial dimension of the mesh vertices
 *
 * \note Assumes that all elements of the mesh have the same number of
 * vertices (as does the MeshData class)
 */
template<int D>
class MeshWrapper
{
private:
   using VertSet = slam::PositionSet<IndexT>;
   using ElemSet = slam::PositionSet<IndexT>;
   using RTStride = slam::policies::RuntimeStride<IndexT>;
   using Card = slam::policies::ConstantCardinality<IndexT, RTStride>;
   using Ind = slam::policies::CArrayIndirection<IndexT, const IndexT>;
   using ElemVertRelation = slam::StaticRelation<IndexT, IndexT, Card, Ind,ElemSet,VertSet>;

public:
   using PointType = primal::Point<double, D>;
   using BBox = primal::BoundingBox<double, D>;

   MeshWrapper() : m_meshData(nullptr) {}

   /*!
    * Constructs MeshWrapper instance from a non-null meshdata pointer
    */
   MeshWrapper(const MeshData* meshData)
      : m_meshData(meshData)
      , m_vertSet(m_meshData->numberOfNodes())
      , m_elemSet(m_meshData->numberOfElements())
   {
      // Generate connectivity relation for elements
      using BuilderType = typename ElemVertRelation::RelationBuilder;

      m_elemVertConnectivity = BuilderType()
          .fromSet( &m_elemSet)
          .toSet( &m_vertSet)
          .begins( typename BuilderType::BeginsSetBuilder()
                   .stride( m_meshData->numberOfNodesPerElement()))
          .indices( typename BuilderType::IndicesSetBuilder()
                    .size( m_elemSet.size() * m_meshData->numberOfNodesPerElement())
                    .data( m_meshData->getConnectivityData() ));
   }

   /*!
    * Gets the vertex at global (host) index \a vId from the mesh
    * \param vId Vertex Id
    * \return A primal Point instance
    */
   PointType getVertex(IndexT vId)
   {
      return PointType::make_point(
            static_cast<double>(m_meshData->getPosition(0)[vId]),
            static_cast<double>(m_meshData->getPosition(1)[vId]),
            (D == 3) ? static_cast<double>(m_meshData->getPosition(2)[vId]) : double() );
   }

   /*!
    * Gets the bounding box for the element with index \a eId
    * \param vId Vertex Id
    * \return A primal BoundingBox instance
    */
   BBox elementBoundingBox(IndexT eId)
   {
      BBox box;

      for(auto vId : m_elemVertConnectivity[eId])
      {
         box.addPoint( getVertex(vId) );
      }

      return box;
   }

   /*! Returns the number of vertices in the mesh */
   int numVerts() const { return m_vertSet.size(); }

   /*! Returns the number of elements in the mesh */
   int numElems() const { return m_elemSet.size(); }

private:
   const MeshData* m_meshData;

   VertSet m_vertSet;
   ElemSet m_elemSet;
   ElemVertRelation m_elemVertConnectivity;
};

/*!
 * \brief Helper class to compute the candidate pairs for a coupling scheme
 *
 * A GridSearch indexes the elements from the first mesh of
 * the coupling scheme in a spatial index that requires element bounding boxes.
 * Then, for each of the elements in the second mesh, we find the candidate
 * pairs and add them to the coupling scheme's list of contact pairs.
 *
 * The spatial index is generated in \a generateSpatialIndex()
 * and the search is performed in \a findInterfacePairs()
 *
 * \tparam D The spatial dimension of the coupling scheme mesh vertices.
 */
template<int D>
class GridSearch
{
public:
   using ImplicitGridType = spin::ImplicitGrid<D, axom::SEQ_EXEC, int>;
   using SpacePoint = typename ImplicitGridType::SpacePoint;
   using SpaceVec = typename ImplicitGridType::SpaceVec;
   using SpatialBoundingBox = typename ImplicitGridType::SpatialBoundingBox;

   /*!
    * Constructs a GridSearch instance over CouplingScheme \a couplingScheme
    * \pre couplingScheme is not null
    */
   GridSearch(CouplingScheme* couplingScheme)
      : m_couplingScheme(couplingScheme)
   {
      MeshManager & meshManager = MeshManager::getInstance();

      IndexT mesh_id1 = m_couplingScheme->getMeshId1();
      MeshData const & meshData1 = *meshManager.at(mesh_id1);
      m_meshWrapper1 = MeshWrapper<D>(&meshData1);

      IndexT mesh_id2 = m_couplingScheme->getMeshId2();
      MeshData const & meshData2 = *meshManager.at(mesh_id2);
      m_meshWrapper2 = MeshWrapper<D>(&meshData2);

      m_couplingScheme->getInterfacePairs()->clear();
   }

   /*!
    * Constructs spatial index over elements of coupling scheme's first mesh
    */
   void generateSpatialIndex()
   {
      // TODO does this tolerance need to scale with the mesh?
      constexpr double bboxTolerance = 1e-6;

      // Find the bounding boxes of the elements in the first mesh
      // Store them in an array for efficient reuse
      m_gridBBox.clear();
      m_meshBBoxes1.reserve(m_meshWrapper1.numElems());
      for(int i=0; i< m_meshWrapper1.numElems(); ++i)
      {
         m_meshBBoxes1.emplace_back(m_meshWrapper1.elementBoundingBox(i));
      }

      // Find an appropriate resolution for the spatial index grid
      //
      // (Note KW) This implementation is a bit ad-hoc
      // * Inflate bounding boxes by 33% of longest dimension
      //   to avoid zero-width dimensions
      // * Find the average extents (range) of the boxes
      //   Assumption is that elements are roughly the same size
      // * Grid resolution for each dimension is overall box width
      //   divided by half the average element width
      SpaceVec ranges;
      for(int i=0; i< m_meshWrapper1.numElems(); ++i)
      {
         auto& bbox = m_meshBBoxes1[i];
         inflateBBox(bbox);

         ranges += bbox.range();

         // build up overall bounding box along the way
         m_gridBBox.addBox( bbox );
      }

      // inflate grid box slightly so elem bounding boxes are not on grid bdry
      m_gridBBox.scale(1.0 + bboxTolerance);

      ranges /= m_meshWrapper1.numElems();

      // Compute grid resolution from average bbox size
      typename ImplicitGridType::GridCell resolution;
      SpaceVec bboxRange = m_gridBBox.range();
      const RealT scaleFac = 0.5; // TODO is this mesh dependent?
      for(int i=0; i < D; ++i)
      {
         resolution[i] = static_cast<IndexT>(
               std::ceil( scaleFac * bboxRange[i] / ranges[i] ));
      }

      // Next, initialize the ImplicitGrid
      m_grid.initialize( m_gridBBox, &resolution, m_meshWrapper1.numElems());

      // Finally, insert the elements
      for(int i=0; i< m_meshWrapper1.numElems(); ++i)
      {
         m_grid.insert(m_meshBBoxes1[i], i);
      }

      // Output some info for debugging
      if(true)
      {
         SLIC_DEBUG("Implicit Grid info: "
             << "\n Mesh 1 bounding box (inflated): " << m_gridBBox
             << "\n Avg range: " << ranges
             << "\n Computed resolution: " << resolution );

         SpatialBoundingBox bbox2;
         for(int i=0; i< m_meshWrapper2.numElems(); ++i)
         {
            bbox2.addBox( m_meshWrapper2.elementBoundingBox(i) );
         }

         SLIC_DEBUG( "Mesh 2 bounding box is: " << bbox2 );
      }

   } // end generateSpatialIndex()

   /*!
    * Use the spatial index to find candidates in first mesh for each
    * element in second mesh of coupling scheme.
    */
   void findInterfacePairs()
   {
      using BitsetType = typename ImplicitGridType::BitsetType;

      // Extract some mesh metadata from coupling scheme / mesh manageer
      MeshManager & meshManager = MeshManager::getInstance();

      IndexT mesh_id1 = m_couplingScheme->getMeshId1();
      MeshData const & meshData1 = *meshManager.at(mesh_id1);
      int cellType1 = static_cast<int>(meshData1.getElementType());

      IndexT mesh_id2 = m_couplingScheme->getMeshId2();
      MeshData const & meshData2 = *meshManager.at(mesh_id2);
      int cellType2 = static_cast<int>(meshData2.getElementType());

      InterfacePairs* contactPairs = m_couplingScheme->getInterfacePairs();


      // Find matches in first mesh (with index 'fromIdx')
      // with candidate elements in second mesh (with index 'toIdx')
      int k = 0;
      for(int toIdx=0; toIdx< m_meshWrapper2.numElems(); ++toIdx)
      {
         SpatialBoundingBox bbox = m_meshWrapper2.elementBoundingBox(toIdx);
         inflateBBox(bbox);

         // Query the mesh
         auto candidateBits = m_grid.getCandidates( bbox );

         // Add candidates
         for(IndexT fromIdx = candidateBits.find_first() ;
             fromIdx != BitsetType::npos ;
             fromIdx = candidateBits.find_next( fromIdx) )
         {
            // if mesh_id1 = mesh_id2, then check to make sure fromIdx < toIdx
            // so we don't double count
            if ( (mesh_id1 == mesh_id2) && (fromIdx < toIdx) )
            {
               continue;
            }

            // TODO: Add extra filter by bbox

            InterfacePair pair( mesh_id1, cellType1, fromIdx,
                                mesh_id2, cellType2, toIdx );

            // perform initial geometry or validity checks to identify initially valid face-pairs
            bool isContactCandidate = geomFilter( pair, m_couplingScheme->getContactMode() );

            // add interface pair for initially valid candidate face-pairs
            if (isContactCandidate)
            {
               pair.pairId = k;
               contactPairs->addInterfacePair( pair );
               ++k;
            }
         }
      }

   } // end findInterfacePairs()


private:
   /*!
    * Expands bounding box by 33% of longest dimension's range
    */
   void inflateBBox(SpatialBoundingBox& bbox)
   {
      constexpr double sc = 1./3.;

      int d = bbox.getLongestDimension();
      const double expansionFac =  sc * bbox.range()[d];
      bbox.expand(expansionFac);
   }

private:

   CouplingScheme* m_couplingScheme;
   MeshWrapper<D> m_meshWrapper1;
   MeshWrapper<D> m_meshWrapper2;

   ImplicitGridType m_grid;
   SpatialBoundingBox m_gridBBox;
   ArrayT<SpatialBoundingBox> m_meshBBoxes1;

};

/*!
 * Compute all pairs of elements in the two meshes of the CouplingScheme
 *
 * \note Assumes the two meshes are different
 */
 void generateCartesianProductPairs(CouplingScheme* cs)
 {
    MeshManager & meshManager = MeshManager::getInstance();

    IndexT mesh_id1 = cs->getMeshId1();
    MeshData const & meshData1 = *meshManager.at(mesh_id1);
    int mesh1NumElems = meshData1.numberOfElements();

    IndexT mesh_id2 = cs->getMeshId2();
    MeshData const & meshData2 = *meshManager.at(mesh_id2);
    int mesh2NumElems = meshData2.numberOfElements();

    int numPairs = mesh1NumElems * mesh2NumElems;

    InterfacePairs* contactPairs = cs->getInterfacePairs();
    contactPairs->clear();
    contactPairs->reserve( numPairs );

    int cellType1 = static_cast<int>(meshData1.getElementType());
    int cellType2 = static_cast<int>(meshData2.getElementType());

    int k = 0;
    for(int fromIdx = 0; fromIdx < mesh1NumElems; ++fromIdx)
    {
       // set starting index for inner loop
       int startIdx = (mesh_id1 == mesh_id2) ? fromIdx : 0;

       for(int toIdx = startIdx; toIdx < mesh2NumElems; ++toIdx)
       {
          InterfacePair pair( mesh_id1, cellType1, fromIdx,
                              mesh_id2, cellType2, toIdx );
          //
          // perform initial geometry or validity checks to identify initially valid face-pairs
          bool isContactCandidate = geomFilter( pair, cs->getContactMode() );

          // add interface pair for initially valid candidate face-pairs
          if (isContactCandidate)
          {
             pair.pairId = k;
             contactPairs->addInterfacePair( pair );
             ++k;
          }
       }
    }

    // print current number of pairs vs. expected for all contact cases other than auto. 
    // the cartesian product will erroneously count a face against itself or faces that share 
    // nodes for auto-contact
    if (cs->getContactCase() != AUTO)
    {
       SLIC_DEBUG("Coupling scheme has " << contactPairs->getNumPairs()
             << " pairs." << " Expected " << numPairs
             << " = " << mesh1NumElems << " * " << mesh2NumElems << ".");
    }

 }



InterfacePairFinder::InterfacePairFinder(CouplingScheme* cs)
   : m_couplingScheme(cs)
   , m_gridSearch2D(nullptr)
   , m_gridSearch3D(nullptr)
{
   SLIC_ASSERT_MSG(m_couplingScheme != nullptr,
         "Coupling scheme was invalid (null pointer)");
}

InterfacePairFinder::~InterfacePairFinder()
{
   if( m_gridSearch2D != nullptr)
   {
      delete m_gridSearch2D;
      m_gridSearch2D = nullptr;
   }

   if( m_gridSearch3D != nullptr)
   {
      delete m_gridSearch3D;
      m_gridSearch3D = nullptr;
   }

}

void InterfacePairFinder::initialize()
{
   const int dim = m_couplingScheme->spatialDimension();

   switch(m_couplingScheme->getBinningMethod() )
   {
   case BINNING_CARTESIAN_PRODUCT:
      // no-op
      break;
   case BINNING_GRID:
      // The spatial grid is templated on the dimension
      switch( dim )
      {
      case 2:
         m_gridSearch2D = new GridSearch<2>(m_couplingScheme);
         m_gridSearch2D->generateSpatialIndex();
         break;
      case 3:
         m_gridSearch3D = new GridSearch<3>(m_couplingScheme);
         m_gridSearch3D->generateSpatialIndex();
         break;
      default:
         SLIC_ERROR("Invalid dimension: " << dim );
         break;
      }
      break;
   default:
      SLIC_ERROR("Unsupported binning method");
      break;
   }
}

void InterfacePairFinder::findInterfacePairs()
{
   const int dim = m_couplingScheme->spatialDimension();

   switch(m_couplingScheme->getBinningMethod() )
   {
   case BINNING_CARTESIAN_PRODUCT:
      generateCartesianProductPairs( m_couplingScheme );
      break;
   case BINNING_GRID:
      // The spatial grid is templated on the dimension
      switch( dim )
      {
      case 2:
         m_gridSearch2D->findInterfacePairs();
         break;
      case 3:
         m_gridSearch3D->findInterfacePairs();
         break;
      default:
         SLIC_ERROR("Invalid dimension: " << dim );
         break;
      }
      break;
   default:
      SLIC_ERROR("Unsupported binning method");
      break;
   }

    // set boolean on coupling scheme object indicating
    // that binning has occurred
    m_couplingScheme->setBinned(true);
}


} // end namespace tribol

