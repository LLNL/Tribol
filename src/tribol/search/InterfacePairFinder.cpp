/*
 *******************************************************************************
 * Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 *******************************************************************************
 */

#include "tribol/search/InterfacePairFinder.hpp"

#include "tribol/common/ExecModel.hpp"
#include "tribol/mesh/CouplingScheme.hpp"
#include "tribol/mesh/MeshManager.hpp"
#include "tribol/mesh/MeshData.hpp"
#include "tribol/mesh/InterfacePairs.hpp"
#include "tribol/utils/Algorithm.hpp"
#include "tribol/utils/Math.hpp"

#include "axom/slic.hpp"
#include "axom/slam.hpp"
#include "axom/primal.hpp"
#include "axom/spin.hpp"

#include "umpire/ResourceManager.hpp"
#include "umpire/TypedAllocator.hpp"
#include "umpire/strategy/DynamicPoolList.hpp"
#include <RAJA/pattern/atomic.hpp>
#include <axom/core/Array.hpp>

// Define some namespace aliases to help with axom usage
namespace slam = axom::slam;
namespace primal = axom::primal;
namespace spin = axom::spin;


namespace tribol
{

/*!
 *  Perform geometry/proximity checks 1-4
 */
TRIBOL_HOST_DEVICE bool geomFilter( IndexT element1, IndexT element2,
                                    const MeshData::Viewer& mesh1, const MeshData::Viewer& mesh2,
                                    ContactMode const mode )
{
  /// CHECK #1: Check to make sure the two face ids are not the same 
  ///           and the two mesh ids are not the same.
  if ((mesh1.meshId() == mesh2.meshId()) && (element1 == element2))
  {
    return false;
  }

  int dim = mesh1.spatialDimension();

  /// CHECK #2: Check to make sure faces don't share a common 
  ///           node for the case where meshId1 = meshId2. 
  ///           We want to preclude two adjacent faces from interacting.
  if (mesh1.meshId() == mesh2.meshId()) 
  {
    for (IndexT i{0}; i < mesh1.numberOfNodesPerElement(); ++i)
    {
      int node1 = mesh1.getGlobalNodeId(element1, i);
      for (IndexT j{0}; j < mesh2.numberOfNodesPerElement(); ++j)
      {
        int node2 = mesh2.getGlobalNodeId(element2, j);
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
      nrmlCheck += mesh1.getElementNormals()[d][element1]
        * mesh2.getElementNormals()[d][element2];
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
      RealT r1 = mesh1.getFaceRadii()[ element1 ];
      RealT r2 = mesh2.getFaceRadii()[ element2 ];

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
      RealT distX = mesh2.getElementCentroids()[0][ element2 ] - mesh1.getElementCentroids()[0][ element1 ];
      RealT distY = mesh2.getElementCentroids()[1][ element2 ] - mesh1.getElementCentroids()[1][ element1 ];
      RealT distZ = mesh2.getElementCentroids()[2][ element2 ] - mesh1.getElementCentroids()[2][ element1 ];
      
      RealT distMag = magnitude(distX, distY, distZ );

      if (distMag >= (distMax)) {
         return false;
      } 
   } // end of dim == 3
   else if (dim == 2)
   {
      // get 1/2 edge length off the mesh data
      RealT e1 = 0.5 * mesh1.getElementAreas()[ element1 ];
      RealT e2 = 0.5 * mesh2.getElementAreas()[ element2 ];

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
      RealT distX = mesh2.getElementCentroids()[0][ element2 ] - mesh1.getElementCentroids()[0][ element1 ];
      RealT distY = mesh2.getElementCentroids()[1][ element2 ] - mesh1.getElementCentroids()[1][ element1 ];

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
   using PointT = primal::Point<RealT, D>;
   using VectorT = primal::Vector<RealT, D>;
   using BBox = primal::BoundingBox<RealT, D>;

   MeshWrapper() : m_meshData(nullptr) {}

   /*!
    * Constructs MeshWrapper instance from a non-null meshdata pointer
    */
   MeshWrapper(const MeshData::Viewer* meshData)
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
   PointT getVertex(IndexT vId)
   {
      return PointT::make_point(
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
      return m_meshData->getFaceRadii()[eId];
   }

   /*!
    * Gets the normal vector for the element with index \a eId
    * \param eId element Id
    * \return A primal Vector instance
    */
   VectorT getFaceNormal(IndexT eId)
   {
      return VectorT::make_vector(
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
   const MeshData::Viewer* m_meshData;

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
      : m_coupling_scheme(couplingScheme)
   {
   }


   void initialize() override 
   {
   }


   void findInterfacePairs() override
   {
      auto& mesh1 = m_coupling_scheme->getMesh1();
      IndexT mesh1NumElems = mesh1.numberOfElements();

      auto& mesh2 = m_coupling_scheme->getMesh2();
      IndexT mesh2NumElems = mesh2.numberOfElements();

      // Reserve memory for boolean array indicating which pairs
      // are in contact
      int numPairs = mesh1NumElems * mesh2NumElems;
      ArrayT<bool> contactArray(numPairs, numPairs, m_coupling_scheme->getAllocatorId());
      bool* inContact = contactArray.data();

      // Allocate memory for a counter
      ArrayT<int> countArray(1, 1, m_coupling_scheme->getAllocatorId());
      int* pCount = countArray.data();

      auto cellType1 = mesh1.getElementType();
      auto cellType2 = mesh2.getElementType();
      ContactMode cmode = m_coupling_scheme->getContactMode();

      IndexT meshId1 = m_coupling_scheme->getMeshId1();
      IndexT meshId2 = m_coupling_scheme->getMeshId2();

      try
      {
         forAllExec(m_coupling_scheme->getExecutionMode(), numPairs,
         [=] TRIBOL_HOST_DEVICE (IndexT i)
         {
            IndexT fromIdx = i / mesh1NumElems;
            IndexT toIdx = i % mesh2NumElems;
            inContact[i] = geomFilter( fromIdx, toIdx, 
                                       mesh1, mesh2,
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

      InterfacePairs* contactPairs = m_coupling_scheme->getInterfacePairs();
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
   CouplingScheme* m_coupling_scheme;
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
      : m_coupling_scheme( couplingScheme )
      , m_mesh_wrapper1( &m_coupling_scheme->getMesh1() )
      , m_mesh_wrapper2( &m_coupling_scheme->getMesh2() )
   {}

   /*!
    * Constructs spatial index over elements of coupling scheme's first mesh
    */
   void initialize() override
   {
      // TODO does this tolerance need to scale with the mesh?
      const RealT bboxTolerance = 1e-6;

      m_coupling_scheme->getInterfacePairs()->clear();

      // Find the bounding boxes of the elements in the first mesh
      // Store them in an array for efficient reuse
      m_gridBBox.clear();
      m_meshBBoxes1.reserve(m_mesh_wrapper1.numElems());
      for(int i=0; i< m_mesh_wrapper1.numElems(); ++i)
      {
         m_meshBBoxes1.emplace_back(m_mesh_wrapper1.elementBoundingBox(i));
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
      for(int i=0; i< m_mesh_wrapper1.numElems(); ++i)
      {
         auto& bbox = m_meshBBoxes1[i];
         inflateBBox(bbox);

         ranges += bbox.range();

         // build up overall bounding box along the way
         m_gridBBox.addBox( bbox );
      }

      // inflate grid box slightly so elem bounding boxes are not on grid bdry
      m_gridBBox.scale(1 + bboxTolerance);

      ranges /= m_mesh_wrapper1.numElems();

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
      m_grid.initialize( m_gridBBox, &resolution, m_mesh_wrapper1.numElems());

      // Finally, insert the elements
      for(int i=0; i< m_mesh_wrapper1.numElems(); ++i)
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
         for(int i=0; i< m_mesh_wrapper2.numElems(); ++i)
         {
            bbox2.addBox( m_mesh_wrapper2.elementBoundingBox(i) );
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
      auto& mesh1 = m_coupling_scheme->getMesh1();
      auto cellType1 = mesh1.getElementType();
      
      auto& mesh2 = m_coupling_scheme->getMesh2();
      auto cellType2 = mesh2.getElementType();

      InterfacePairs* contactPairs = m_coupling_scheme->getInterfacePairs();

      IndexT meshId1 = m_coupling_scheme->getMeshId1();
      IndexT meshId2 = m_coupling_scheme->getMeshId2();

      // Find matches in first mesh (with index 'fromIdx')
      // with candidate elements in second mesh (with index 'toIdx')
      int k = 0;
      for(int toIdx=0; toIdx< m_mesh_wrapper2.numElems(); ++toIdx)
      {
         SpatialBoundingBox bbox = m_mesh_wrapper2.elementBoundingBox(toIdx);
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
            if ( (mesh1.meshId() == mesh2.meshId()) && (fromIdx < toIdx) )
            {
               continue;
            }

            // TODO: Add extra filter by bbox

            // Preliminary geometry/proximity checks, SRW
            bool contact = geomFilter( fromIdx, toIdx,
                                       mesh1, mesh2,
                                       m_coupling_scheme->getContactMode() );

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

   CouplingScheme* m_coupling_scheme;
   MeshWrapper<D> m_mesh_wrapper1;
   MeshWrapper<D> m_mesh_wrapper2;

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
  using BVHT = axom::spin::BVH<D, ExecSpace, RealT>;
  using BoxT = typename BVHT::BoxType;
  using PointT = primal::Point<RealT, D>;
  using RayT = primal::Ray<RealT, D>;
  using VectorT = primal::Vector<RealT, D>;
  using AtomicPolicy = typename axom::execution_space<ExecSpace>::atomic_policy;

  /*!
  * Constructs a BvhSearch instance over CouplingScheme \a couplingScheme
  * \pre couplingScheme is not null
  */
  BvhSearch(CouplingScheme* coupling_scheme)
    : m_coupling_scheme( coupling_scheme )
    , m_mesh_wrapper1( &m_coupling_scheme->getMesh1() )
    , m_mesh_wrapper2( &m_coupling_scheme->getMesh2() )
    , m_boxes1( axom::ArrayOptions::Uninitialized{},
                m_mesh_wrapper1.numElems(),
                m_mesh_wrapper1.numElems(),
                m_coupling_scheme->getAllocatorId() )
    , m_boxes2( axom::ArrayOptions::Uninitialized{},
                m_mesh_wrapper2.numElems(),
                m_mesh_wrapper2.numElems(),
                m_coupling_scheme->getAllocatorId() )
    , m_candidates( axom::ArrayOptions::Uninitialized{},
                    0,
                    0,
                    m_coupling_scheme->getAllocatorId() )
    , m_offsets( axom::ArrayOptions::Uninitialized{},
                 m_mesh_wrapper2.numElems(),
                 m_mesh_wrapper2.numElems(),
                 m_coupling_scheme->getAllocatorId() )
    , m_counts( axom::ArrayOptions::Uninitialized{},
                m_mesh_wrapper2.numElems(),
                m_mesh_wrapper2.numElems(),
                m_coupling_scheme->getAllocatorId() )
  {}

  /*!
  * Allocate and fill bounding box arrays for each of the two meshes
  */
  void initialize() override
  {
    buildMeshBBoxes(m_boxes1, m_coupling_scheme->getMesh1());
    buildMeshBBoxes(m_boxes2, m_coupling_scheme->getMesh2());
  } // end initialize()
   

  /*!
  * Use the BVH to find candidates in first mesh for each
  * element in second mesh of coupling scheme.
  */
  void findInterfacePairs() override
  {
    // Build the BVH
    BVHT bvh;
    bvh.setAllocatorID(m_coupling_scheme->getAllocatorId());
    bvh.initialize(m_boxes1.view(), m_boxes1.size());

    // Search for intersecting bounding boxes
    bvh.findBoundingBoxes(m_offsets.view(),
                          m_counts.view(),
                          m_candidates,
                          m_mesh_wrapper2.numElems(),
                          m_boxes2.view());

    // Apply geom filter to check if intersecting bounding boxes are in contact
    // Change candidate value to -1 if geom filter checks are failed
    auto counts_view = m_counts.view();
    auto offsets_view = m_offsets.view();
    auto candidates_view = m_candidates.view();
    ArrayT<IndexT> filtered_candidates_data(1, 1, m_coupling_scheme->getAllocatorId());
    auto filtered_candidates = filtered_candidates_data.view();
    auto& mesh1 = m_coupling_scheme->getMesh1();
    auto& mesh2 = m_coupling_scheme->getMesh2();
    auto cmode = m_coupling_scheme->getContactMode();
    forAllExec(m_coupling_scheme->getExecutionMode(), m_candidates.size(),
      [=] TRIBOL_HOST_DEVICE (IndexT i) 
      {
        auto mesh1_elem = algorithm::binarySearch(offsets_view, counts_view, i);
        auto mesh2_elem = candidates_view[i];
        if (geomFilter(mesh1_elem, mesh2_elem, mesh1, mesh2, cmode))
        {
          RAJA::atomicInc<AtomicPolicy>(filtered_candidates.data());
        }
        else
        {
          candidates_view[i] = -1;
        }
      }
    );

    // Add filtered pairs to interface pairs
    ArrayT<IndexT> filtered_candidates_host( filtered_candidates_data, 
                                             getResourceAllocatorID(MemorySpace::Host) );
    m_coupling_scheme->getInterfacePairs()->resize(filtered_candidates_host[0]);
    filtered_candidates_host[0] = 0;
    filtered_candidates_data = filtered_candidates_host;

    m_coupling_scheme->getInterfacePairs()->setMeshId(1, m_coupling_scheme->getMeshId1());
    m_coupling_scheme->getInterfacePairs()->setMeshId(2, m_coupling_scheme->getMeshId2());
    m_coupling_scheme->getInterfacePairs()->setPairType(
      1, m_coupling_scheme->getMesh1().getElementType());
    m_coupling_scheme->getInterfacePairs()->setPairType(
      2, m_coupling_scheme->getMesh2().getElementType());

    filtered_candidates = filtered_candidates_data.view();
    auto idx1_view = m_coupling_scheme->getInterfacePairs()->getPairIndex1Array();
    auto idx2_view = m_coupling_scheme->getInterfacePairs()->getPairIndex2Array();
    auto contact_view = m_coupling_scheme->getInterfacePairs()->getContactArray();
    forAllExec(m_coupling_scheme->getExecutionMode(), m_candidates.size(),
      [=] TRIBOL_HOST_DEVICE (IndexT i)
      {
        // Filtering removed this case
        if (candidates_view[i] == -1)
        {
          return;
        }
        
        auto mesh1_elem = algorithm::binarySearch(offsets_view, counts_view, i);
        auto mesh2_elem = candidates_view[i];

        // get unique index for the array
        auto idx = RAJA::atomicInc<AtomicPolicy>(filtered_candidates.data());

        idx1_view[idx] = mesh1_elem;
        idx2_view[idx] = mesh2_elem;
        contact_view[idx] = true;
      }
    );
   } // end findInterfacePairs()

  void buildMeshBBoxes(ArrayT<BoxT>& boxes, const MeshData::Viewer& mesh)
  {
    auto boxes1_view = boxes.view();
    forAllExec(m_coupling_scheme->getExecutionMode(), mesh.numberOfElements(),
      [=] TRIBOL_HOST_DEVICE (IndexT i) {
        BoxT box;
        auto num_nodes_per_elem = mesh.numberOfNodesPerElement();
        for(IndexT j{0}; j < num_nodes_per_elem; ++j)
        {
          IndexT node_id = mesh.getGlobalNodeId(i, j);
          RealT pos[3];
          pos[0] = mesh.getPosition()[0][node_id];
          pos[1] = mesh.getPosition()[1][node_id];
          pos[2] = mesh.getPosition()[2][node_id];  // unused if D==2
          box.addPoint( PointT(pos) );
        }
        // Expand the bounding box in the face normal direction
        RealT vnorm[3];
        mesh.getFaceNormal(i, vnorm);
        VectorT faceNormal(vnorm);
        RealT faceRadius = mesh.getFaceRadii()[i];
        expandBBoxNormal(box, faceNormal, faceRadius);
        boxes1_view[i] = std::move(box);
      }
    );
  }

private:
  /*!
  * Expands bounding box by projecting the face normal by a distance 
  * equal to the effective face radius
  */
  TRIBOL_HOST_DEVICE void expandBBoxNormal(BoxT& bbox, 
                                          const VectorT& faceNormal, 
                                          const RealT faceRadius)
  {
    PointT p0 = bbox.getCentroid();
    RayT outwardRay(p0, faceNormal);
    VectorT inwardNormal(faceNormal);
    inwardNormal *= -1.0;  // this operation is available on device
    RayT inwardRay(p0, inwardNormal);
    PointT pout = outwardRay.at(faceRadius);
    PointT pin = inwardRay.at(faceRadius);
    bbox.addPoint(pout);
    bbox.addPoint(pin);
  }

   /*!
    * Isotropically expands bounding box by the effective face radius.
    */
   TRIBOL_HOST_DEVICE void inflateBBox(BoxT& bbox, 
                                       const RealT faceRadius)
   {
      bbox.expand(faceRadius);
   }

  CouplingScheme* m_coupling_scheme;
  MeshWrapper<D> m_mesh_wrapper1;
  MeshWrapper<D> m_mesh_wrapper2;
  ArrayT<BoxT> m_boxes1;
  ArrayT<BoxT> m_boxes2;
  ArrayT<IndexT> m_candidates;
  ArrayT<IndexT> m_offsets;
  ArrayT<IndexT> m_counts;
};  // End of BvhSearch class definition

///////////////////////////////////////////////////////////////////////////////

InterfacePairFinder::InterfacePairFinder(CouplingScheme* cs)
   : m_coupling_scheme(cs)
{
   SLIC_ASSERT_MSG(cs != nullptr, "Coupling scheme was invalid (null pointer)");
   const int dim = m_coupling_scheme->spatialDimension();
   axom::slic::flushStreams();
   m_search = nullptr;
   switch(cs->getBinningMethod() )
   {
   case BINNING_CARTESIAN_PRODUCT:
      switch( dim )
      {
      case 2:
         m_search = new CartesianProduct<2>(m_coupling_scheme);
         break;
      case 3:
         m_search = new CartesianProduct<3>(m_coupling_scheme);
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
         m_search = new GridSearch<2>(m_coupling_scheme);
         break;
      case 3:
         m_search = new GridSearch<3>(m_coupling_scheme);
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
               m_search = new BvhSearch<2, axom::SEQ_EXEC>(m_coupling_scheme);
               break;
            #ifdef TRIBOL_USE_OPENMP
            case(ExecutionMode::OpenMP):  // This causes compiler to hang
               //SLIC_ERROR("Unsupported execution mode: " << parameters.exec_mode );
               m_search = new BvhSearch<2, axom::OMP_EXEC>(m_coupling_scheme);
               break;
            #endif
            #ifdef TRIBOL_USE_CUDA
            case(ExecutionMode::Cuda):
               m_search = new BvhSearch<2, axom::CUDA_EXEC<TRIBOL_BLOCK_SIZE>>(m_coupling_scheme);
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
               m_search = new BvhSearch<3, axom::SEQ_EXEC>(m_coupling_scheme);
               break;
            #ifdef TRIBOL_USE_OPENMP
            case(ExecutionMode::OpenMP): // This causes compiler to hang
               //SLIC_ERROR("Unsupported execution mode: " << parameters.exec_mode );
               m_search = new BvhSearch<3, axom::OMP_EXEC>(m_coupling_scheme);
               break;
            #endif
            #ifdef TRIBOL_USE_CUDA
            case(ExecutionMode::Cuda):
               m_search = new BvhSearch<3, axom::CUDA_EXEC<TRIBOL_BLOCK_SIZE>>(m_coupling_scheme);
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
   m_coupling_scheme->setBinned(true);
}


} // end namespace tribol

