/*
 *******************************************************************************
 * Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 *******************************************************************************
 */

#include "tribol/search/InterfacePairFinder.hpp"

#include "tribol/types.hpp"
#include "tribol/mesh/CouplingScheme.hpp"
#include "tribol/mesh/MeshManager.hpp"
#include "tribol/mesh/MeshData.hpp"
#include "tribol/mesh/InterfacePairs.hpp"
#include "tribol/common/Parameters.hpp"
#include "tribol/common/loop_exec.hpp"
#include "tribol/utils/Math.hpp"

#include "axom/slic.hpp"
#include "axom/slam.hpp"
#include "axom/primal.hpp"
#include "axom/spin.hpp"
#include "axom/core/execution/execution_space.hpp" 

#include "umpire/ResourceManager.hpp"
#include "umpire/TypedAllocator.hpp"

// Define some namespace aliases to help with axom usage
namespace slam = axom::slam;
namespace primal = axom::primal;
namespace spin = axom::spin;


namespace tribol
{

/*!
 *  Perform geometry/proximity checks 1-4
 */
TRIBOL_HOST_DEVICE bool geomFilter( const integer pairIndex1, const integer pairIndex2,
                                    const integer meshId1, const integer meshId2, 
                                    const MeshData* const pMesh1, const MeshData* const pMesh2,
                                    ContactMode const mode )
{
   /// CHECK #1: Check to make sure the two face ids are not the same 
   ///           and the two mesh ids are not the same.
   if ((meshId1 == meshId2) && (pairIndex1 == pairIndex2))
   {
      return false;
   }

   int dim = ( pMesh1->m_elementType == tribol::EDGE ) ? 2 : 3;

   /// CHECK #2: Check to make sure faces don't share a common 
   ///           node for the case where meshId1 = meshId2. 
   ///           We want to preclude two adjacent faces from interacting.
   if (meshId1 == meshId2) 
   {
      for (int i=0; i<pMesh1->m_numCellNodes; ++i)
      {
         int node1 = pMesh1->getFaceNodeId(pairIndex1, i);
         for (int j=0; j<pMesh2->m_numCellNodes; ++j)
         {
            int node2 = pMesh2->getFaceNodeId(pairIndex2, j);
            if (node1 == node2)
            {
              return false;
            }
         }
      }
   }

   /// CHECK #3: Check that face normals are opposing up to some tolerance.
   ///           This uses a hard coded normal tolerance for this check.
   real nrmlTol = -0.173648177; // taken as cos(100) between face pair

   real m_nZ1, m_nZ2;
   if ( dim == 3 )
   {
      m_nZ1 = pMesh1->m_nZ[ pairIndex1 ];
      m_nZ2 = pMesh2->m_nZ[ pairIndex2 ];
   }
   else
   {
      m_nZ1 = 0.;
      m_nZ2 = 0.;
   }

   real nrmlCheck = pMesh1->m_nX[pairIndex1] * pMesh2->m_nX[pairIndex2] + 
                    pMesh1->m_nY[pairIndex1] * pMesh2->m_nY[pairIndex2] + 
                    m_nZ1 * m_nZ2;

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
   real offset_tol = 0.05;
   if (dim == 3)
   {
      real r1 = pMesh1->m_faceRadius[ pairIndex1 ];
      real r2 = pMesh2->m_faceRadius[ pairIndex2 ];

      // set maximum offset of face centroids for inclusion
      real distMax = r1 + r2; // default is sum of face radii

      // check if the contact mode is conforming, in which case the 
      // faces are supposed to be aligned
      if (mode == SURFACE_TO_SURFACE_CONFORMING)
      { 
         // use 5% of max face radius for conforming case as 
         // tolerance on face offsets
         distMax *= offset_tol; 
      }

      // compute the distance between the two face centroids
      real distX = pMesh2->m_cX[ pairIndex2 ] - pMesh1->m_cX[ pairIndex1 ];
      real distY = pMesh2->m_cY[ pairIndex2 ] - pMesh1->m_cY[ pairIndex1 ];
      real distZ = pMesh2->m_cZ[ pairIndex2 ] - pMesh1->m_cZ[ pairIndex1 ];
      
      real distMag = magnitude(distX, distY, distZ );

      if (distMag >= (distMax)) {
         return false;
      } 
   } // end of dim == 3
   else if (dim == 2)
   {
      // get 1/2 edge length off the mesh data
      real e1 = 0.5 * pMesh1->m_area[ pairIndex1 ];
      real e2 = 0.5 * pMesh2->m_area[ pairIndex2 ];

      // set maximum offset of edge centroids for inclusion
      real distMax = e1 + e2; // default is sum of 1/2 edge lengths

      // check if the contact mode is conforming, in which case the 
      // edges are supposed to be aligned
      if (mode == SURFACE_TO_SURFACE_CONFORMING)
      { 
         // use 5% of max face radius for conforming case as 
         // tolerance on face offsets
         distMax *= offset_tol; 
      }

      // compute the distance between the two edge centroids
      real distX = pMesh2->m_cX[ pairIndex2 ] - pMesh1->m_cX[ pairIndex1 ];
      real distY = pMesh2->m_cY[ pairIndex2 ] - pMesh1->m_cY[ pairIndex1 ];

      real distMag = magnitude(distX, distY);

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
   using VertSet = slam::PositionSet<IndexType>;
   using ElemSet = slam::PositionSet<IndexType>;
   using RTStride = slam::policies::RuntimeStride<IndexType>;
   using Card = slam::policies::ConstantCardinality<IndexType, RTStride>;
   using Ind = slam::policies::ArrayIndirection<IndexType, const IndexType>;
   using ElemVertRelation = slam::StaticRelation<IndexType, IndexType, Card, Ind,ElemSet,VertSet>;

public:
   using PointType = primal::Point<real, D>;
   using VectorType = primal::Vector<real, D>;
   using BBox = primal::BoundingBox<real, D>;

   MeshWrapper() : m_meshData(nullptr) {}

   /*!
    * Constructs MeshWrapper instance from a non-null meshdata pointer
    */
   MeshWrapper(const MeshData* meshData)
      : m_meshData(meshData)

      , m_vertSet(m_meshData->m_lengthNodalData)
      , m_elemSet(m_meshData->m_numCells)
   {
      // Generate connectivity relation for elements
      using BuilderType = typename ElemVertRelation::RelationBuilder;

      m_elemVertConnectivity = BuilderType()
          .fromSet( &m_elemSet)
          .toSet( &m_vertSet)
          .begins( typename BuilderType::BeginsSetBuilder()
                   .stride( m_meshData->m_numCellNodes))
          .indices( typename BuilderType::IndicesSetBuilder()
                    .size( m_elemSet.size() * m_meshData->m_numCellNodes)
                    .data( m_meshData->m_connectivity ));
   }

   /*!
    * Gets the vertex at global (host) index \a vId from the mesh
    * \param vId Vertex Id
    * \return A primal Point instance
    */
   PointType getVertex(IndexType vId)
   {
      return PointType::make_point(
            m_meshData->m_positionX[vId],
            m_meshData->m_positionY[vId],
            (D == 3) ? m_meshData->m_positionZ[vId] : real() );
   }

   /*!
    * Returns the area for the element with index \a eId
    * \param eId element Id
    * \return A primal Vector instance
    */
   real getFaceArea(IndexType eId)
   {
      return m_meshData->m_area[eId];
   }

   /*!
    * Returns the effective radius of the element with index \a eId
    * \param eId element Id
    * \return face area
    */
   real getFaceRadius(IndexType eId)
   {
      return m_meshData->m_faceRadius[eId];
   }

   /*!
    * Gets the normal vector for the element with index \a eId
    * \param eId element Id
    * \return A primal Vector instance
    */
   VectorType getFaceNormal(IndexType eId)
   {
      return VectorType::make_vector(
            m_meshData->m_nX[eId],
            m_meshData->m_nY[eId],
            (D == 3) ? m_meshData->m_nZ[eId] : real() );
   }

   /*!
    * Gets the bounding box for the element with index \a eId
    * \param eId element Id
    * \return A primal BoundingBox instance
    */
   BBox elementBoundingBox(IndexType eId)
   {
      BBox box;

      for(auto vId : m_elemVertConnectivity[eId])
      {
         box.addPoint( getVertex(vId) );
      }
      return box;
   }


   /*! Returns the number of vertices in the mesh */
   integer numVerts() const { return m_vertSet.size(); }

   /*! Returns the number of elements in the mesh */
   integer numElems() const { return m_elemSet.size(); }

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
      MeshManager& meshManager = MeshManager::getInstance();

      integer meshId1 = m_couplingScheme->getMeshId1();
      const MeshData* pMesh1 = &meshManager.GetMeshInstance(meshId1);
      integer mesh1NumElems = pMesh1->m_numCells;

      integer meshId2 = m_couplingScheme->getMeshId2();
      const MeshData* pMesh2 = &meshManager.GetMeshInstance(meshId2); // TODO: make const
      integer mesh2NumElems = pMesh2->m_numCells;

      // Reserve memory for boolean array indicating which pairs
      // are in contact
      int numPairs = mesh1NumElems * mesh2NumElems;
      Array1D<bool> contactArray;  
      contactArray.resize(numPairs);
      bool* inContact = contactArray.data();

      // Allocate memory for a counter
      Array1D<int> countArray;  
      countArray.resize(1);
      countArray[0] = 0;
      int* pCount = countArray.data();

      int cellType1 = static_cast<integer>(pMesh1->m_elementType);
      int cellType2 = static_cast<integer>(pMesh2->m_elementType);
      ContactMode cmode = m_couplingScheme->getContactMode();

      try
      {
         TRIBOL_FORALL(k, numPairs, 
         {
            IndexType fromIdx = k / mesh1NumElems;
            IndexType toIdx = k % mesh2NumElems;
            inContact[k] = geomFilter( fromIdx, toIdx, 
                                       meshId1, meshId2,
                                       pMesh1, pMesh2,
                                       cmode );
            RAJA::atomicAdd< RAJA::auto_atomic >(pCount, static_cast<int>(inContact[k]));
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
         SCOPED_RANGE("addInterfacePairs", 4);
         for (int k=0; k<numPairs; k++) 
         {
            if (inContact[k])
            {
               IndexType fromIdx = k / mesh1NumElems;
               IndexType toIdx = k % mesh2NumElems;
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
   using ImplicitGridType = spin::ImplicitGrid<D,int>;
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

      integer meshId1 = m_couplingScheme->getMeshId1();
      const MeshData& meshData1 = meshManager.GetMeshInstance(meshId1);
      m_meshWrapper1 = MeshWrapper<D>(&meshData1);

      integer meshId2 = m_couplingScheme->getMeshId2();
      const MeshData& meshData2 = meshManager.GetMeshInstance(meshId2);
      m_meshWrapper2 = MeshWrapper<D>(&meshData2);
   }

   /*!
    * Constructs spatial index over elements of coupling scheme's first mesh
    */
   void initialize() override
   {
      // TODO does this tolerance need to scale with the mesh?
      const real bboxTolerance = 1e-6;

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
      const real scaleFac = 0.5; // TODO is this mesh dependent?
      for(int i=0; i < D; ++i)
      {
         resolution[i] = static_cast<IndexType>(
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
      MeshManager & meshManager = MeshManager::getInstance();

      integer meshId1 = m_couplingScheme->getMeshId1();
      const MeshData& meshData1 = meshManager.GetMeshInstance(meshId1);
      const MeshData* pMesh1 = &meshData1;
      int cellType1 = static_cast<integer>(meshData1.m_elementType);

      integer meshId2 = m_couplingScheme->getMeshId2();
      const MeshData& meshData2 = meshManager.GetMeshInstance(meshId2);
      const MeshData* pMesh2 = &meshData2;
      int cellType2 = static_cast<integer>(meshData2.m_elementType);

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
         for(IndexType fromIdx = candidateBits.find_first() ;
             fromIdx != BitsetType::npos ;
             fromIdx = candidateBits.find_next( fromIdx) )
         {
            // if meshId1 = meshId2, then check to make sure fromIdx < toIdx 
            // so we don't double count
            if ( (meshId1 == meshId2) && (fromIdx < toIdx) )
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
      const real expansionFac =  sc * bbox.range()[d];
      bbox.expand(expansionFac);
   }

   CouplingScheme* m_couplingScheme;
   MeshWrapper<D> m_meshWrapper1;
   MeshWrapper<D> m_meshWrapper2;

   ImplicitGridType m_grid;
   SpatialBoundingBox m_gridBBox;
   Array1D<SpatialBoundingBox> m_meshBBoxes1;

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
   using PointType = primal::Point<real, D>;
   using RayType = primal::Ray<real, D>;
   using VectorType = primal::Vector<real, D>;
   using BoxType = primal::BoundingBox<real, D>;
   using BVHType = spin::BVH<D, ExecSpace, real>;

   /*!
    * Constructs a BvhSearch instance over CouplingScheme \a couplingScheme
    * \pre couplingScheme is not null
    */
   BvhSearch(CouplingScheme* couplingScheme)
      : m_couplingScheme(couplingScheme)
      , m_allocatorId(0)
      , m_boxes1(nullptr)
      , m_boxes2(nullptr)
      , m_offsets(nullptr)
      , m_counts(nullptr)
      , m_candidates(nullptr)
      , m_contact(nullptr)
      , m_idx(nullptr)
   {
      MeshManager & meshManager = MeshManager::getInstance();

      integer meshId1 = m_couplingScheme->getMeshId1();
      MeshData const & meshData1 = meshManager.GetMeshInstance(meshId1);
      m_mesh1 = &meshData1;
      m_meshWrapper1 = MeshWrapper<D>(&meshData1);

      integer meshId2 = m_couplingScheme->getMeshId2();
      MeshData const & meshData2 = meshManager.GetMeshInstance(meshId2);
      m_mesh2 = &meshData2;
      m_meshWrapper2 = MeshWrapper<D>(&meshData2);
   }

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
         allocator.deallocate(m_candidates);
         allocator.deallocate(m_contact);
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
         rm.makeAllocator<umpire::strategy::DynamicPool>(
         allocator.getName() + "_POOL",
         allocator,
         POOL_SIZE);
      m_allocatorId = pool_allocator.getId();
      umpire::TypedAllocator<BoxType> boxAllocator{pool_allocator};

      // Find the bounding boxes of all elements in the first mesh
      // and store them in m_boxes1.  
      const int N1 = m_mesh1->m_numCells;
      const int NUM_CELL_NODES1 = m_mesh1->m_numCellNodes;
      m_boxes1 = boxAllocator.allocate(N1);

      #ifdef __NVCC__   // Workaround for compiler issues
      TRIBOL_FORALL(i, N1, 
      #else
      using EXEC_POL = typename axom::execution_space<ExecSpace>::loop_policy;
      RAJA::forall<EXEC_POL>(RAJA::RangeSegment(0, N1), [=] (int i)
      #endif
      {
         BoxType box;
         for(int j=0; j<NUM_CELL_NODES1; ++j)
         {
            int nodeIndex = (NUM_CELL_NODES1 * i) + j;
            int nodeId = m_mesh1->m_connectivity[ nodeIndex ];
            real pos[3];
            pos[0] = m_mesh1->m_positionX[nodeId];
            pos[1] = m_mesh1->m_positionY[nodeId];
            pos[2] = m_mesh1->m_positionZ[nodeId];  // unused if D==2
            box.addPoint( PointType(pos) );
         }
         // Expand the bounding box in the face normal direction
         real vnorm[3];
         vnorm[0] = m_mesh1->m_nX[i];
         vnorm[1] = m_mesh1->m_nY[i];
         vnorm[2] = m_mesh1->m_nZ[i];   // unused if D==2
         VectorType faceNormal = VectorType(vnorm);
         real faceRadius = m_mesh1->m_faceRadius[i];
         expandBBoxNormal(box, faceNormal, faceRadius);
         m_boxes1[i] = box;
      });

      // Find the bounding boxes of all elements in the second mesh
      // and store them in m_boxes2.  
      const int N2 = m_mesh2->m_numCells;
      const int NUM_CELL_NODES2 = m_mesh2->m_numCellNodes;
      m_boxes2 = boxAllocator.allocate(N2);
      #ifdef __NVCC__   // Workaround for compiler issues
      TRIBOL_FORALL(i, N2, 
      #else
      RAJA::forall<EXEC_POL>(RAJA::RangeSegment(0, N2), [=] (int i)
      #endif
      {
         BoxType box;
         for(int j=0; j<NUM_CELL_NODES2; ++j)
         {
            int nodeIndex = (NUM_CELL_NODES2 * i) + j;
            int nodeId = m_mesh2->m_connectivity[ nodeIndex ];
            real pos[3];
            pos[0] = m_mesh2->m_positionX[nodeId];
            pos[1] = m_mesh2->m_positionY[nodeId];
            pos[2] = m_mesh2->m_positionZ[nodeId];  // unused if D==2
            box.addPoint( PointType(pos) );
         }
         // Use scale() since expand() is host-only
         real faceRadius = m_mesh2->m_faceRadius[i];
         real vlen = box.range().norm();
         real factor = 1.0 + faceRadius / vlen;
         box.scale(factor);
         m_boxes2[i] = box;
      });

      // Allocate offset and count arrays. The BVH query results will 
      // be stored here in findInterfacePairs().
      umpire::TypedAllocator<IndexType> idxAllocator{pool_allocator};
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
      MeshManager & meshManager = MeshManager::getInstance();

      integer meshId1 = m_couplingScheme->getMeshId1();
      const MeshData* pMesh1 = &meshManager.GetMeshInstance(meshId1);
      int cellType1 = static_cast<integer>(pMesh1->m_elementType);

      integer meshId2 = m_couplingScheme->getMeshId2();
      const MeshData* pMesh2 = &meshManager.GetMeshInstance(meshId2);
      int cellType2 = static_cast<integer>(pMesh2->m_elementType);

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
         SCOPED_RANGE("bvh.initialize", 8);
         auto status = bvh.initialize(m_boxes1, N1);
         SLIC_ERROR_IF((status != spin::BVH_BUILD_OK), "BVH initialization failed");
      }
 
      // Output some diagnostic information
      BoxType bounds = bvh.getBounds();
      VectorType pmin = bounds.getMin();
      VectorType pmax = bounds.getMax();
      SLIC_INFO("BVH info: "
            << "\n Mesh 1 expanded root bounding box axis minima: " << pmin
            << "\n Mesh 1 expanded root bounding box axis maxima: " << pmax );

      // Query the BVH to find intersections between elements of mesh 1
      // and mesh2.
      {
         SCOPED_RANGE("findBoundingBoxes", 9);
         bvh.findBoundingBoxes(m_offsets, m_counts, m_candidates, N2, m_boxes2);
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
      #ifdef __NVCC__   // Workaround for compiler issues
      TRIBOL_FORALL(i, N2, 
      #else
      RAJA::forall<EXEC_POL>(RAJA::RangeSegment(0, N2), [=] (int i)
      #endif
      {
         candTotal += static_cast<int>(m_counts[i]); 
      });
      const int ncand = static_cast<int>(candTotal.get());
      SLIC_INFO(" Found " << ncand << " candidate pairs");
      if (ncand < 1)
      {
         return;
      }
      SLIC_ERROR_IF((m_candidates == nullptr), "Empty candidates array");
      auto& rm = umpire::ResourceManager::getInstance();
      umpire::Allocator allocator = rm.getAllocator(m_allocatorId);
      umpire::TypedAllocator<IndexType> idxAllocator{allocator};
      m_contact = idxAllocator.allocate(ncand);
      m_idx = idxAllocator.allocate(ncand);

      // Find candidate intersecting bounding boxes in first mesh 
      // (with index 'fromIdx') for elements in second mesh (with index 'toIdx')
      const ContactMode cmode = m_couplingScheme->getContactMode();
      {
         SCOPED_RANGE("geomFilter", 10);
         // FIXME: This kernel is susceptible to thread divergence
         #ifdef __NVCC__   // Workaround for compiler issues
         TRIBOL_FORALL(toIdx, N2, 
         #else
         RAJA::forall<EXEC_POL>(RAJA::RangeSegment(0, N2), [=] (int toIdx)
         #endif
         {
            for (IndexType k=0; k<m_counts[toIdx]; k++) 
            {
               IndexType pairIdx = m_offsets[toIdx] + k;
               IndexType fromIdx = m_candidates[pairIdx];
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
      RAJA::exclusive_scan<EXEC_POL>(m_contact, 
                                     m_contact+ncand,
                                     m_idx, 
                                     RAJA::operators::plus<IndexType>{});

      const IndexType contact_total = m_idx[ncand-1] + 1;
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
      integer* const pidx1 = contactPairs->getPairIndex1Array();
      integer* const pidx2 = contactPairs->getPairIndex2Array();
      bool* const contact = contactPairs->getContactArray();
      #ifdef __NVCC__   // Workaround for compiler issues
      TRIBOL_FORALL(toIdx, N2, 
      #else
      RAJA::forall<EXEC_POL>(RAJA::RangeSegment(0, N2), [=] (int toIdx)
      #endif
      {
         for (IndexType k=0; k<m_counts[toIdx]; k++) 
         {
               IndexType pairIdx = m_offsets[toIdx] + k;
               IndexType fromIdx = m_candidates[pairIdx];
               IndexType idx = m_idx[pairIdx];
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
                                            const real faceRadius)
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
                                       const real faceRadius)
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
   IndexType* m_offsets;
   IndexType* m_counts;
   IndexType* m_candidates;
   IndexType* m_contact;
   IndexType* m_idx;
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
         switch(parameters.exec_mode)
         {
            case(LoopExecMode::SEQUENTIAL):
               m_search = new BvhSearch<2, axom::SEQ_EXEC>(m_couplingScheme);
               break;
            #ifdef AXOM_USE_OPENMP
            case(LoopExecMode::OPENMP_PARALLEL):  // This causes compiler to hang
               //SLIC_ERROR("Unsupported execution mode: " << parameters.exec_mode );
               m_search = new BvhSearch<2, axom::OMP_EXEC>(m_couplingScheme);
               break;
            #endif
            #ifdef AXOM_USE_CUDA
            case(LoopExecMode::CUDA_PARALLEL):
               m_search = new BvhSearch<2, axom::CUDA_EXEC<CUDA_BLOCK_SIZE>>(m_couplingScheme);
               break;
            #endif
            default:
               SLIC_ERROR("Invalid execution mode: " << parameters.exec_mode );
               break;
         }
         break;
      case 3:
         switch(parameters.exec_mode)
         {
            case(LoopExecMode::SEQUENTIAL):
               m_search = new BvhSearch<3, axom::SEQ_EXEC>(m_couplingScheme);
               break;
            #ifdef AXOM_USE_OPENMP
            case(LoopExecMode::OPENMP_PARALLEL): // This causes compiler to hang
               //SLIC_ERROR("Unsupported execution mode: " << parameters.exec_mode );
               m_search = new BvhSearch<3, axom::OMP_EXEC>(m_couplingScheme);
               break;
            #endif
            #ifdef AXOM_USE_CUDA
            case(LoopExecMode::CUDA_PARALLEL):
               m_search = new BvhSearch<3, axom::CUDA_EXEC<CUDA_BLOCK_SIZE>>(m_couplingScheme);
               break;
            #endif
            default:
               SLIC_ERROR("Invalid execution mode: " << parameters.exec_mode );
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
   SCOPED_RANGE("initialize", 6);
   m_search->initialize();
}

void InterfacePairFinder::findInterfacePairs()
{
   SLIC_INFO("Searching for interface pairs");
   SCOPED_RANGE("findInterfacePairs", 7);
   m_search->findInterfacePairs();
   // set boolean on coupling scheme object indicating 
   // that binning has occurred
   m_couplingScheme->setBinned(true);
}


} // end namespace tribol

