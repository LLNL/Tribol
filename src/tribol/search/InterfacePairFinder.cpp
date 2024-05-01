/*
 *******************************************************************************
 * Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 *******************************************************************************
 */

#include "tribol/search/InterfacePairFinder.hpp"

#include "tribol/mesh/CouplingScheme.hpp"
#include "tribol/mesh/MeshManager.hpp"
#include "tribol/mesh/MeshData.hpp"
#include "tribol/mesh/InterfacePairs.hpp"
#include "tribol/utils/Math.hpp"

#include "axom/slic.hpp"
#include "axom/slam.hpp"
#include "axom/primal.hpp"
#include "axom/spin.hpp"

#include "umpire/ResourceManager.hpp"
#include "umpire/TypedAllocator.hpp"
#include "umpire/strategy/DynamicPoolList.hpp"

// Define some namespace aliases to help with axom usage
namespace slam = axom::slam;
namespace primal = axom::primal;
namespace spin = axom::spin;


namespace tribol
{

/*!
 *  Perform geometry/proximity checks 1-4
 */
TRIBOL_HOST_DEVICE bool geomFilter( const IndexT pairIndex1, const IndexT pairIndex2,
                                    const IndexT meshId1, const IndexT meshId2, 
                                    const MeshData* const pMesh1, const MeshData* const pMesh2,
                                    ContactMode const mode )
{
   /// CHECK #1: Check to make sure the two face ids are not the same 
   ///           and the two mesh ids are not the same.
   if ((meshId1 == meshId2) && (pairIndex1 == pairIndex2))
   {
      return false;
   }

   int dim = pMesh1->dimension();

   /// CHECK #2: Check to make sure faces don't share a common 
   ///           node for the case where meshId1 = meshId2. 
   ///           We want to preclude two adjacent faces from interacting.
   if (meshId1 == meshId2) 
   {
      for (int i=0; i<pMesh1->numberOfNodesPerElement(); ++i)
      {
         int node1 = pMesh1->getGlobalNodeId(pairIndex1, i);
         for (int j=0; j<pMesh2->numberOfNodesPerElement(); ++j)
         {
            int node2 = pMesh2->getGlobalNodeId(pairIndex2, j);
            if (node1 == node2)
            {
              return false;
            }
         }
      }
   }

   /// CHECK #3: Check that face normals are opposing up to some tolerance.
   ///           This uses a hard coded normal tolerance for this check.
   RealT nrmlTol = -0.173648177; // taken as cos(100) between face pair
   
   RealT nrmlCheck = 0.0;
   for (int d{0}; d < dim; ++d)
   {
      nrmlCheck += pMesh1->getElementNormals()[d][pairIndex1]
        * pMesh2->getElementNormals()[d][pairIndex2];
   }

   // check normal projection against tolerance
   if (nrmlCheck > nrmlTol) {
      return false;
   }

   /// CHECK #4 (3D): Perform radius check, which involves seeing if 
   ///                the distance between the two face vertex averaged
   ///                centroid is less than the sum of the two face radii.
   ///                The face radii are taken to be the magnitude of the 
   ///                longest vector from that face's vertex averaged 
   ///                centroid to one its nodes.
   RealT offset_tol = 0.05;
   if (dim == 3)
   {
      RealT r1 = pMesh1->getFaceRadius()[ pairIndex1 ];
      RealT r2 = pMesh2->getFaceRadius()[ pairIndex2 ];

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
      RealT distX = pMesh2->getElementCentroids()[0][ pairIndex2 ] - pMesh1->getElementCentroids()[0][ pairIndex1 ];
      RealT distY = pMesh2->getElementCentroids()[1][ pairIndex2 ] - pMesh1->getElementCentroids()[1][ pairIndex1 ];
      RealT distZ = pMesh2->getElementCentroids()[2][ pairIndex2 ] - pMesh1->getElementCentroids()[2][ pairIndex1 ];
      
      RealT distMag = magnitude(distX, distY, distZ );

      if (distMag >= (distMax)) {
         return false;
      } 
   } // end of dim == 3
   else if (dim == 2)
   {
      // get 1/2 edge length off the mesh data
      RealT e1 = 0.5 * pMesh1->getElementAreas()[ pairIndex1 ];
      RealT e2 = 0.5 * pMesh2->getElementAreas()[ pairIndex2 ];

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
      RealT distX = pMesh2->getElementCentroids()[0][ pairIndex2 ] - pMesh1->getElementCentroids()[0][ pairIndex1 ];
      RealT distY = pMesh2->getElementCentroids()[1][ pairIndex2 ] - pMesh1->getElementCentroids()[1][ pairIndex1 ];

      RealT distMag = magnitude(distX, distY);

      if (distMag >= (distMax)) 
      {
         return false;
      }
   } // end of dim == 2

   // if we made it here we passed all checks
   return true;

} // end geomFilter()


//------------------------------------------------------------------------------


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
   using Ind = slam::policies::ArrayIndirection<IndexT, IndexT>;
   using ElemVertRelation = slam::StaticRelation<IndexT, IndexT, Card, Ind,ElemSet,VertSet>;

public:
   using PointType = primal::Point<RealT, D>;
   using VectorType = primal::Vector<RealT, D>;
   using BBox = primal::BoundingBox<RealT, D>;

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
      axom::Array<IndexT> conn;

      m_elemVertConnectivity = BuilderType()
          .fromSet( &m_elemSet)
          .toSet( &m_vertSet)
          .begins( typename BuilderType::BeginsSetBuilder()
                   .stride( m_meshData->numberOfNodesPerElement()))
          .indices( typename BuilderType::IndicesSetBuilder()
                    .size( m_elemSet.size() * m_meshData->numberOfNodesPerElement())
                    .data( &conn ));
   }

   /*!
    * Gets the vertex at global (host) index \a vId from the mesh
    * \param vId Vertex Id
    * \return A primal Point instance
    */
   PointType getVertex(IndexT vId)
   {
      return PointType::make_point(
            m_meshData->getPosition()[0][vId],
            m_meshData->getPosition()[1][vId],
            (D == 3) ? m_meshData->getPosition()[2][vId] : RealT() );
   }

   /*!
    * Returns the area for the element with index \a eId
    * \param eId element Id
    * \return A primal Vector instance
    */
   RealT getFaceArea(IndexT eId)
   {
      return m_meshData->getElementAreas()[eId];
   }

   /*!
    * Returns the effective radius of the element with index \a eId
    * \param eId element Id
    * \return face area
    */
   RealT getFaceRadius(IndexT eId)
   {
      return m_meshData->getFaceRadius()[eId];
   }

   /*!
    * Gets the normal vector for the element with index \a eId
    * \param eId element Id
    * \return A primal Vector instance
    */
   VectorType getFaceNormal(IndexT eId)
   {
      return VectorType::make_vector(
            m_meshData->getElementNormals()[0][eId],
            m_meshData->getElementNormals()[1][eId],
            (D == 3) ? m_meshData->getElementNormals()[2][eId] : RealT() );
   }

   /*!
    * Gets the bounding box for the element with index \a eId
    * \param eId element Id
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
   IndexT numVerts() const { return m_vertSet.size(); }

   /*! Returns the number of elements in the mesh */
   IndexT numElems() const { return m_elemSet.size(); }

private:
   const MeshData* m_meshData;

   VertSet m_vertSet;
   ElemSet m_elemSet;
   ElemVertRelation m_elemVertConnectivity;
};

/*!
 * \brief Base class to compute the candidate pairs for a coupling scheme
 *
 * \a initialize() must be called prior to \a findInterfacePairs()
 *
 */
class SearchBase
{
   public:
   //using SpacePoint = typename ImplicitGridType::SpacePoint;
   //using SpaceVec = typename ImplicitGridType::SpaceVec;
   //using SpatialBoundingBox = typename ImplicitGridType::SpatialBoundingBox;
   SearchBase() {};
   virtual ~SearchBase() {};
   /*!
    * Prepares the object for spatial searches
    */
   virtual void initialize() = 0;

   /*!
    * Find candidates in first mesh for each element in second mesh of coupling scheme.
    */
   virtual void findInterfacePairs() = 0;
};

///////////////////////////////////////////////////////////////////////////////

/*!
 * \brief Helper class to compute the candidate pairs for a coupling scheme
 *
 * A CartesianProduct search combines each element from the first mesh of
 * the coupling scheme with each element in the second mesh. A geometry filter
 * is then applied to each resulting element pair.  This is the slowest of all
 * pair-finding methods since ALL possible element pairs are considered, i.e.,
 * this is an exhaustive search.
 *
 * \tparam D The spatial dimension of the coupling scheme mesh vertices.
 */
template<int D>
class CartesianProduct : public SearchBase
{
public:
   /*!
    * Constructs a CartesianProduct instance over CouplingScheme \a couplingScheme
    * \pre couplingScheme is not null
    */
   CartesianProduct(CouplingScheme* couplingScheme)
      : m_couplingScheme(couplingScheme)
   {
   }


   void initialize() override 
   {
   }


   void findInterfacePairs() override
   {
      const MeshData* pMesh1 = m_couplingScheme->getMesh1();
      IndexT mesh1NumElems = pMesh1->numberOfElements();

      const MeshData* pMesh2 = m_couplingScheme->getMesh2();
      IndexT mesh2NumElems = pMesh2->numberOfElements();

      // Reserve memory for boolean array indicating which pairs
      // are in contact
      int numPairs = mesh1NumElems * mesh2NumElems;
      ArrayT<bool> contactArray(numPairs, numPairs, m_couplingScheme->getAllocatorId());
      bool* inContact = contactArray.data();

      // Allocate memory for a counter
      ArrayT<int> countArray(1, 1, m_couplingScheme->getAllocatorId());
      int* pCount = countArray.data();

      int cellType1 = static_cast<IndexT>(pMesh1->getElementType());
      int cellType2 = static_cast<IndexT>(pMesh2->getElementType());
      ContactMode cmode = m_couplingScheme->getContactMode();

      IndexT meshId1 = m_couplingScheme->getMeshId1();
      IndexT meshId2 = m_couplingScheme->getMeshId2();

      try
      {
         forAllExec(m_couplingScheme->getExecutionMode(), numPairs,
         [mesh1NumElems, mesh2NumElems, inContact, meshId1, meshId2, pMesh1, pMesh2, cmode, pCount] TRIBOL_HOST_DEVICE (const IndexT i)
         {
            IndexT fromIdx = i / mesh1NumElems;
            IndexT toIdx = i % mesh2NumElems;
            inContact[i] = geomFilter( fromIdx, toIdx, 
                                       meshId1, meshId2,
                                       pMesh1, pMesh2,
                                       cmode );
            RAJA::atomicAdd< RAJA::auto_atomic >(pCount, static_cast<int>(inContact[i]));
         });
      }  // End of profiling block
      catch(const std::exception& e)
      {
         std::cerr << e.what() << std::endl;
      }
      catch(...)
      {
         std::cerr << "unknown exception encountered during findInterfacePairs/geomFilter" << std::endl;
      }
      
      SLIC_INFO("Found " << countArray[0] << " pairs in contact" );

      InterfacePairs* contactPairs = m_couplingScheme->getInterfacePairs();
      contactPairs->clear();
      contactPairs->reserve(countArray[0]);
      
      int idx = 0;
      {
         for (int k=0; k<numPairs; k++) 
         {
            if (inContact[k])
            {
               IndexT fromIdx = k / mesh1NumElems;
               IndexT toIdx = k % mesh2NumElems;
               InterfacePair pair( meshId1, cellType1, fromIdx,
                                   meshId2, cellType2, toIdx, 
                                   true, idx );
               // SLIC_INFO("Interface pair " << idx << " = " << toIdx << ", " << fromIdx);  // Debug only
               contactPairs->addInterfacePair( pair );
               idx++;
            }
         }
      }
      SLIC_INFO("Added " << idx << " pairs to contact list" );  // Debug only

      SLIC_INFO("Coupling scheme has " << contactPairs->getNumPairs()
            << " pairs out of a maximum possible of " << numPairs
            << " = " << mesh1NumElems << " * " << mesh2NumElems << ".");
   }
private:
   CouplingScheme* m_couplingScheme;
};  // End of CartesianProduct definition

///////////////////////////////////////////////////////////////////////////////

/*!
 * \brief Implicit Grid helper class to compute the candidate pairs for a coupling scheme
 *
 * A GridSearch indexes the elements from the first mesh of
 * the coupling scheme in a spatial index that requires element bounding boxes.
 * Then, for each of the elements in the second mesh, we find the candidate
 * pairs and add them to the coupling scheme's list of contact pairs.
 *
 * The spatial index is generated in \a initialize()
 * and the search is performed in \a findInterfacePairs()
 *
 * \tparam D The spatial dimension of the coupling scheme mesh vertices.
 */
template<int D>
class GridSearch : public SearchBase
{
public:
   using ImplicitGridType = spin::ImplicitGrid<D,axom::SEQ_EXEC,int>;
   using SpacePoint = typename ImplicitGridType::SpacePoint;
   using SpaceVec = typename ImplicitGridType::SpaceVec;
   using SpatialBoundingBox = typename ImplicitGridType::SpatialBoundingBox;

   /*!
    * Constructs a GridSearch instance over CouplingScheme \a couplingScheme
    * \pre couplingScheme is not null
    */
   GridSearch( CouplingScheme* couplingScheme )
      : m_couplingScheme( couplingScheme )
      , m_meshWrapper1( m_couplingScheme->getMesh1() )
      , m_meshWrapper2( m_couplingScheme->getMesh2() )
   {}

   /*!
    * Constructs spatial index over elements of coupling scheme's first mesh
    */
   void initialize() override
   {
      // TODO does this tolerance need to scale with the mesh?
      const RealT bboxTolerance = 1e-6;

      m_couplingScheme->getInterfacePairs()->clear();

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
      m_gridBBox.scale(1 + bboxTolerance);

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
         SLIC_INFO("Implicit Grid info: "
             << "\n Mesh 1 bounding box (inflated): " << m_gridBBox
             << "\n Avg range: " << ranges
             << "\n Computed resolution: " << resolution );

         SpatialBoundingBox bbox2;
         for(int i=0; i< m_meshWrapper2.numElems(); ++i)
         {
            bbox2.addBox( m_meshWrapper2.elementBoundingBox(i) );
         }

         SLIC_INFO( "Mesh 2 bounding box is: " << bbox2 );
      }
   }; // end initialize()

   /*!
    * Use the spatial index to find candidates in first mesh for each
    * element in second mesh of coupling scheme.
    */
   void findInterfacePairs() override
   {
      using BitsetType = typename ImplicitGridType::BitsetType;

      // Extract some mesh metadata from coupling scheme / mesh manageer
      const MeshData* pMesh1 = m_couplingScheme->getMesh1();
      int cellType1 = static_cast<IndexT>(pMesh1->getElementType());
      
      const MeshData* pMesh2 = m_couplingScheme->getMesh2();
      int cellType2 = static_cast<IndexT>(pMesh2->getElementType());

      InterfacePairs* contactPairs = m_couplingScheme->getInterfacePairs();

      IndexT meshId1 = m_couplingScheme->getMeshId1();
      IndexT meshId2 = m_couplingScheme->getMeshId2();

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
            // if meshId1 = meshId2, then check to make sure fromIdx < toIdx 
            // so we don't double count
            if ( (pMesh1 == pMesh2) && (fromIdx < toIdx) )
            {
               continue;
            }

            // TODO: Add extra filter by bbox

            // Preliminary geometry/proximity checks, SRW
            bool contact = geomFilter( fromIdx, toIdx, 
                                       meshId1, meshId2,
                                       pMesh1, pMesh2,
                                       m_couplingScheme->getContactMode() );

            if (contact)
            {
               InterfacePair pair( meshId1, cellType1, fromIdx,
                                 meshId2, cellType2, toIdx, false, 
                                 -1 );
               pair.pairId = k;
               contactPairs->addInterfacePair( pair );
               // SLIC_INFO("Interface pair " << k << " = " << toIdx << ", " << fromIdx);  // Debug only
               ++k;
            }
         }
      }  // end of loop over candidates in second mesh

   } // end findInterfacePairs()


private:
   /*!
    * Expands bounding box by 33% of longest dimension's range
    */
   void inflateBBox(SpatialBoundingBox& bbox)
   {
      constexpr double sc = 1./3.;

      int d = bbox.getLongestDimension();
      const RealT expansionFac =  sc * bbox.range()[d];
      bbox.expand(expansionFac);
   }

   CouplingScheme* m_couplingScheme;
   MeshWrapper<D> m_meshWrapper1;
   MeshWrapper<D> m_meshWrapper2;

   ImplicitGridType m_grid;
   SpatialBoundingBox m_gridBBox;
   ArrayT<SpatialBoundingBox> m_meshBBoxes1;

}; // End of GridSearch class definition


///////////////////////////////////////////////////////////////////////////////

/*!
 * \brief BVH helper class to compute the candidate pairs for a coupling scheme
 *
 * A BvhSearch indexes the elements from the first mesh of
 * the coupling scheme in a spatial index that requires element bounding boxes.
 * Then, for each of the elements in the second mesh, we find the candidate
 * pairs and add them to the coupling scheme's list of contact pairs.
 *
 * The spatial index is generated in \a generateSpatialIndex()
 * and the search is performed in \a findInterfacePairs()
 * 
 * \tparam D The spatial dimension of the coupling scheme mesh vertices.
 */
template<int D, class ExecSpace>
class BvhSearch : public SearchBase
{
public:
   using PointType = primal::Point<RealT, D>;
   using RayType = primal::Ray<RealT, D>;
   using VectorType = primal::Vector<RealT, D>;
   using BoxType = primal::BoundingBox<RealT, D>;
   using BVHType = spin::BVH<D, ExecSpace, RealT>;

   /*!
    * Constructs a BvhSearch instance over CouplingScheme \a couplingScheme
    * \pre couplingScheme is not null
    */
   BvhSearch(CouplingScheme* couplingScheme)
      : m_couplingScheme(couplingScheme)
      , m_meshWrapper1(m_couplingScheme->getMesh1())
      , m_meshWrapper2(m_couplingScheme->getMesh2())
      , m_allocatorId(0)
      , m_boxes1(nullptr)
      , m_boxes2(nullptr)
      , m_offsets(nullptr)
      , m_counts(nullptr)
      , m_idx(nullptr)
   {}

   /*!
    * Clean up
    */
   ~BvhSearch() 
   {
      if (m_allocatorId != 0)
      {
         auto& rm = umpire::ResourceManager::getInstance();
         umpire::Allocator allocator = rm.getAllocator(m_allocatorId);
         allocator.deallocate(m_boxes1);
         allocator.deallocate(m_boxes2);
         allocator.deallocate(m_offsets);
         allocator.deallocate(m_counts);
         allocator.deallocate(m_idx);
      }
   }

   /*!
    * Allocate and fill bounding box arrays for each of the two meshes
    */
   void initialize() override
   {
      const int ndim = m_couplingScheme->spatialDimension();
      SLIC_ERROR_IF( (ndim != D), D << "-dimensional BvhSearch created for "<< ndim << "-dimensional coupling scheme");
      m_couplingScheme->getInterfacePairs()->clear();
      constexpr size_t POOL_SIZE = (1024 * 1024 * 1024) + 1;
      auto& rm = umpire::ResourceManager::getInstance();
      auto resource = (axom::execution_space<ExecSpace>::onDevice()) ? umpire::resource::Unified : umpire::resource::Host;
      umpire::Allocator allocator = rm.getAllocator(resource);
      umpire::Allocator pool_allocator =
         rm.makeAllocator<umpire::strategy::DynamicPoolList>(
         allocator.getName() + "_POOL",
         allocator,
         POOL_SIZE);
      m_allocatorId = pool_allocator.getId();
      umpire::TypedAllocator<BoxType> boxAllocator{pool_allocator};

      // Find the bounding boxes of all elements in the first mesh
      // and store them in m_boxes1.  
      const int N1 = m_mesh1->numberOfElements();
      const int NUM_CELL_NODES1 = m_mesh1->numberOfNodesPerElement();
      m_boxes1 = boxAllocator.allocate(N1);

      using EXEC_POL = typename axom::execution_space<ExecSpace>::loop_policy;
      RAJA::forall<EXEC_POL>(RAJA::RangeSegment(0, N1), [=] (int i)
      {
         BoxType box;
         for(int j=0; j<NUM_CELL_NODES1; ++j)
         {
            int nodeIndex = (NUM_CELL_NODES1 * i) + j;
            int nodeId = m_mesh1->getConnectivity().data()[ nodeIndex ];
            RealT pos[3];
            pos[0] = m_mesh1->getPosition()[0][nodeId];
            pos[1] = m_mesh1->getPosition()[1][nodeId];
            pos[2] = m_mesh1->getPosition()[2][nodeId];  // unused if D==2
            box.addPoint( PointType(pos) );
         }
         // Expand the bounding box in the face normal direction
         RealT vnorm[3];
         vnorm[0] = m_mesh1->getElementNormals()[0][i];
         vnorm[1] = m_mesh1->getElementNormals()[1][i];
         vnorm[2] = m_mesh1->getElementNormals()[2][i];   // unused if D==2
         VectorType faceNormal = VectorType(vnorm);
         RealT faceRadius = m_mesh1->getFaceRadius()[i];
         expandBBoxNormal(box, faceNormal, faceRadius);
         m_boxes1[i] = box;
      });

      // Find the bounding boxes of all elements in the second mesh
      // and store them in m_boxes2.  
      const int N2 = m_mesh2->numberOfElements();
      const int NUM_CELL_NODES2 = m_mesh2->numberOfNodesPerElement();
      m_boxes2 = boxAllocator.allocate(N2);
      RAJA::forall<EXEC_POL>(RAJA::RangeSegment(0, N2), [=] (int i)
      {
         BoxType box;
         for(int j=0; j<NUM_CELL_NODES2; ++j)
         {
            int nodeIndex = (NUM_CELL_NODES2 * i) + j;
            int nodeId = m_mesh2->getConnectivity().data()[ nodeIndex ];
            RealT pos[3];
            pos[0] = m_mesh2->getPosition()[0][nodeId];
            pos[1] = m_mesh2->getPosition()[1][nodeId];
            pos[2] = m_mesh2->getPosition()[2][nodeId];  // unused if D==2
            box.addPoint( PointType(pos) );
         }
         // Use scale() since expand() is host-only
         RealT faceRadius = m_mesh2->getFaceRadius()[i];
         RealT vlen = box.range().norm();
         RealT factor = 1.0 + faceRadius / vlen;
         box.scale(factor);
         m_boxes2[i] = box;
      });

      // Allocate offset and count arrays. The BVH query results will 
      // be stored here in findInterfacePairs().
      umpire::TypedAllocator<IndexT> idxAllocator{pool_allocator};
      m_offsets = idxAllocator.allocate(N2);
      m_counts = idxAllocator.allocate(N2);
      rm.memset(m_offsets, 0);      // FIXME: redundant
      rm.memset(m_counts, 0);       // FIXME: redundant
   } // end initialize()
   

   /*!
    * Use the BVH to find candidates in first mesh for each
    * element in second mesh of coupling scheme.
    */
   void findInterfacePairs() override
   {
      SLIC_ERROR_IF((m_allocatorId == 0), "BvhSearch not initialized");

      // Extract some mesh metadata from coupling scheme / mesh manageer
      const MeshData* pMesh1 = m_couplingScheme->getMesh1();
      int cellType1 = static_cast<IndexT>(pMesh1->getElementType());

      const MeshData* pMesh2 = m_couplingScheme->getMesh2();
      int cellType2 = static_cast<IndexT>(pMesh2->getElementType());

      const int N1 = m_meshWrapper1.numElems();
      const int N2 = m_meshWrapper2.numElems();
      // #define TRIBOL_DEBUG_BVH 1
      #ifdef TRIBOL_DEBUG_BVH
      const int N1_OUTPUT = (N1<=16) ? N1 : 16;
      std::cout << "  Mesh 1 bounding boxes" << std::endl;
      for (int i=0; i<N1_OUTPUT; i++)
      {
         std::cout << i << ": " << m_boxes1[i] << std::endl;
      }
      const int N2_OUTPUT = (N2<=16) ? N2 : 16;
      std::cout << std::endl << "  Mesh 2 bounding boxes" << std::endl;
      for (int i=0; i<N2_OUTPUT; i++)
      {
         std::cout << i << ": " << m_boxes2[i] << std::endl;
      }
      #endif

      // Build the BVH
      BVHType bvh;
      bvh.setAllocatorID(m_allocatorId);
      bvh.setScaleFactor(1.01);  // Expand bounding volume dimensions by 1% - probably not needed
      {
         auto status = bvh.initialize(m_boxes1, N1);
         SLIC_ERROR_IF((status != spin::BVH_BUILD_OK), "BVH initialization failed");
      }
 
      // Output some diagnostic information
      BoxType bounds = bvh.getBounds();
      auto pmin = bounds.getMin();
      auto pmax = bounds.getMax();
      SLIC_INFO("BVH info: "
            << "\n Mesh 1 expanded root bounding box axis minima: " << pmin
            << "\n Mesh 1 expanded root bounding box axis maxima: " << pmax );

      // Query the BVH to find intersections between elements of mesh 1
      // and mesh2.
      {
         bvh.findBoundingBoxes(axom::ArrayView<IndexT>(m_offsets, N2), axom::ArrayView<IndexT>(m_counts, N2), m_candidates, N2, m_boxes2);
      }

      #ifdef TRIBOL_DEBUG_BVH
      std::cout << std::endl << "  BVH query results" << std::endl;
      std::cout << "index, counts, candidates" << std::endl;
      for (int i=0; i<N2_OUTPUT; i++)
      {
         std::cout << i << ", " << m_counts[i] << ", " << m_offsets[i];
         for (int j=0; j<m_counts[i]; j++)
         {
            std::cout << ", " << m_candidates[j];
         }
         std::cout << std::endl;
      }
      #endif

      // Get the total number of candidates
      // NOTES: 1. This gives different results if OMP_THREADS>1 on TOSS3 gcc-8.3.1
      //        2. RAJA::forall<EXEC_POL>(RAJA::RangeSegment(0, N2), [=] (int i) gives a compiler error for xl/nvcc
      //  Below is a workaround to compile differently
      //
      using EXEC_POL = typename axom::execution_space<ExecSpace>::loop_policy;
      using REDUCE_POL = typename axom::execution_space<ExecSpace>::reduce_policy;
      RAJA::ReduceSum<REDUCE_POL, int> candTotal(0);
      RAJA::forall<EXEC_POL>(RAJA::RangeSegment(0, N2), [=] (int i)
      {
         candTotal += static_cast<int>(m_counts[i]); 
      });
      const int ncand = static_cast<int>(candTotal.get());
      SLIC_INFO(" Found " << ncand << " candidate pairs");
      if (ncand < 1)
      {
         return;
      }
      SLIC_ERROR_IF((m_candidates.empty()), "Empty candidates array");
      auto& rm = umpire::ResourceManager::getInstance();
      umpire::Allocator allocator = rm.getAllocator(m_allocatorId);
      umpire::TypedAllocator<IndexT> idxAllocator{allocator};
      m_contact.reserve(ncand);
      m_idx = idxAllocator.allocate(ncand);

      IndexT meshId1 = m_couplingScheme->getMeshId1();
      IndexT meshId2 = m_couplingScheme->getMeshId2();

      // Find candidate intersecting bounding boxes in first mesh 
      // (with index 'fromIdx') for elements in second mesh (with index 'toIdx')
      const ContactMode cmode = m_couplingScheme->getContactMode();
      {
         // FIXME: This kernel is susceptible to thread divergence
         RAJA::forall<EXEC_POL>(RAJA::RangeSegment(0, N2), [=] (int toIdx)
         {
            for (IndexT k=0; k<m_counts[toIdx]; k++) 
            {
               IndexT pairIdx = m_offsets[toIdx] + k;
               IndexT fromIdx = m_candidates[pairIdx];
               m_contact[pairIdx] = geomFilter( fromIdx, toIdx, 
                                                meshId1, meshId2,
                                                pMesh1, pMesh2,
                                                cmode );
            }
         });  // end loop over second mesh elements (toIdx)
      }  // end geomFilter profiling range

      /////////////////////////////////////////////////////////////////////////
      // FIXME: The remainder of this function could be factored out
      //        and used by other binning methods to create the list of
      //        contact pairs
      
      // Perform scan to get indexes of contact pairs
      // m_idx is incremented for each true value of m_contact so
      // that the filtered pair index is given by m_idx[pairIdx].
      // Note that m_idx is the cumulative sum of m_contact.
      // 
      // RAJA::exclusive_scan<EXEC_POL>(m_contact, 
      //                                m_contact,
      //                                m_idx);

      const IndexT contact_total = m_idx[ncand-1] + 1;
      SLIC_INFO("Found " << contact_total << " interface pairs in contact" );

      #ifdef TRIBOL_DEBUG_BVH
      const int NCAND_OUTPUT = (ncand<=100) ? ncand : 100;
      std::cout << std::endl << "  Filtered Candidates" << std::endl;
      std::cout << "candidate pair index, contact, index" << std::endl;
      for (int i=0; i<NCAND_OUTPUT; i++)
      {
         std::cout << i << ", " << m_contact[i] << ", " << m_idx[i] << std::endl;
      }
      #endif
      // Create list of contact pairs in the coupling scheme
      // The first and second pairs are recorded in pidx1 and pidx2 
      // if m_contact is true.  Both pidx1 and pidx2 are first resized
      // to length contact_total before filling.
      //
      InterfacePairs* contactPairs = m_couplingScheme->getInterfacePairs();
      contactPairs->setMeshId(meshId1, meshId2);
      contactPairs->setPairType(1, cellType1);
      contactPairs->setPairType(2, cellType2);
      contactPairs->resize(contact_total);
      IndexT* const pidx1 = contactPairs->getPairIndex1Array();
      IndexT* const pidx2 = contactPairs->getPairIndex2Array();
      bool* const contact = contactPairs->getContactArray();
      RAJA::forall<EXEC_POL>(RAJA::RangeSegment(0, N2), [=] (int toIdx)
      {
         for (IndexT k=0; k<m_counts[toIdx]; k++) 
         {
               IndexT pairIdx = m_offsets[toIdx] + k;
               IndexT fromIdx = m_candidates[pairIdx];
               IndexT idx = m_idx[pairIdx];
               if (m_contact[pairIdx])
               {
                  pidx2[idx] = toIdx;
                  pidx1[idx] = fromIdx;
                  contact[idx] = true;
               }
         }
      });  
   } // end findInterfacePairs()


private:
   /*!
    * Expands bounding box by projecting the face normal by a distance 
    * equal to the effective face radius
    */
   TRIBOL_HOST_DEVICE void expandBBoxNormal(BoxType& bbox, 
                                            const VectorType& faceNormal, 
                                            const RealT faceRadius)
   {
      PointType p0 = bbox.getCentroid();
      RayType outwardRay(p0, faceNormal);
      VectorType inwardNormal(faceNormal);
      inwardNormal *= -1.0;  // this operation is available on device
      RayType inwardRay(p0, inwardNormal);
      PointType pout = outwardRay.at(faceRadius);
      PointType pin = inwardRay.at(faceRadius);
      bbox.addPoint(pout);
      bbox.addPoint(pin);
   }

   /*!
    * Isotropically expands bounding box by the effective face radius.
    */
   TRIBOL_HOST_DEVICE void inflateBBox(BoxType& bbox, 
                                       const RealT faceRadius)
   {
      bbox.expand(faceRadius);
   }

   CouplingScheme* m_couplingScheme;
   const MeshData* m_mesh1;
   const MeshData* m_mesh2;
   MeshWrapper<D> m_meshWrapper1;
   MeshWrapper<D> m_meshWrapper2;
   int m_allocatorId;
   BoxType* m_boxes1;
   BoxType* m_boxes2;
   IndexT* m_offsets;
   IndexT* m_counts;
   axom::Array<IndexT> m_candidates;
   axom::Array<IndexT> m_contact;
   IndexT* m_idx;
};  // End of BvhSearch class definition

///////////////////////////////////////////////////////////////////////////////

InterfacePairFinder::InterfacePairFinder(CouplingScheme* cs)
   : m_couplingScheme(cs)
{
   SLIC_ASSERT_MSG(cs != nullptr, "Coupling scheme was invalid (null pointer)");
   constexpr int CUDA_BLOCK_SIZE = 256;
   const int dim = m_couplingScheme->spatialDimension();
   parameters_t& parameters = parameters_t::getInstance();
   axom::slic::flushStreams();
   m_search = nullptr;
   switch(cs->getBinningMethod() )
   {
   case BINNING_CARTESIAN_PRODUCT:
      switch( dim )
      {
      case 2:
         m_search = new CartesianProduct<2>(m_couplingScheme);
         break;
      case 3:
         m_search = new CartesianProduct<3>(m_couplingScheme);
         break;
      default:
         SLIC_ERROR("Invalid dimension: " << dim );
         break;
      } // end of BINNING_CARTESIAN_PRODUCT dimension switch
      break;
   case BINNING_GRID:
      // The spatial grid is templated on the dimension
      switch( dim )
      {
      case 2:
         m_search = new GridSearch<2>(m_couplingScheme);
         break;
      case 3:
         m_search = new GridSearch<3>(m_couplingScheme);
         break;
      default:
         SLIC_ERROR("Invalid dimension: " << dim );
         break;
      } // end of BINNING_GRID dimension switch
      break;
   case BINNING_BVH:
      // The BVH is templated on the dimension and execution space
      switch( dim )
      {
      case 2:
         switch(cs->getExecutionMode())
         {
            case(ExecutionMode::Sequential):
               m_search = new BvhSearch<2, axom::SEQ_EXEC>(m_couplingScheme);
               break;
            #ifdef TRIBOL_USE_OPENMP
            case(ExecutionMode::OpenMP):  // This causes compiler to hang
               //SLIC_ERROR("Unsupported execution mode: " << parameters.exec_mode );
               m_search = new BvhSearch<2, axom::OMP_EXEC>(m_couplingScheme);
               break;
            #endif
            #ifdef AXOM_USE_CUDA
            case(ExecutionMode::Cuda):
               m_search = new BvhSearch<2, axom::CUDA_EXEC<CUDA_BLOCK_SIZE>>(m_couplingScheme);
               break;
            #endif
            default:
               SLIC_ERROR("Invalid execution mode.");
               break;
         }
         break;
      case 3:
         switch(cs->getExecutionMode())
         {
            case(ExecutionMode::Sequential):
               m_search = new BvhSearch<3, axom::SEQ_EXEC>(m_couplingScheme);
               break;
            #ifdef AXOM_USE_OPENMP
            case(ExecutionMode::OpenMP): // This causes compiler to hang
               //SLIC_ERROR("Unsupported execution mode: " << parameters.exec_mode );
               m_search = new BvhSearch<3, axom::OMP_EXEC>(m_couplingScheme);
               break;
            #endif
            #ifdef AXOM_USE_CUDA
            case(ExecutionMode::Cuda):
               m_search = new BvhSearch<3, axom::CUDA_EXEC<CUDA_BLOCK_SIZE>>(m_couplingScheme);
               break;
            #endif
            default:
               SLIC_ERROR("Invalid execution mode.");
               break;
         }
         break;
      default:
         SLIC_ERROR("Invalid dimension: " << dim );
         break;
      } // end of BINNING_BVH dimension switch
      break;
   default:
      SLIC_ERROR("Invalid binning method: " << cs->getBinningMethod() );
      break;
   }  // end of binning method switch
}

InterfacePairFinder::~InterfacePairFinder()
{
   if( m_search != nullptr)
   {
      delete m_search;
   }
}

void InterfacePairFinder::initialize()
{
   SLIC_ASSERT(m_search != nullptr);
   m_search->initialize();
}

void InterfacePairFinder::findInterfacePairs()
{
   SLIC_INFO("Searching for interface pairs");
   m_search->findInterfacePairs();
   // set boolean on coupling scheme object indicating 
   // that binning has occurred
   m_couplingScheme->setBinned(true);
}


} // end namespace tribol

