// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Tribol Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (MIT)

#ifndef SRC_UTILS_TESTUTILS_HPP_
#define SRC_UTILS_TESTUTILS_HPP_

#include "tribol/types.hpp"
#include "tribol/common/Parameters.hpp"

//MFEM includes
#include "mfem.hpp"

#include "axom/config.hpp"


namespace tribol
{

/// RAII struct to initialize and finalize MPI for tests and examples
struct SimpleMPIWrapper
{
   SimpleMPIWrapper(int argc, char* argv[])
   {
   #ifdef TRIBOL_USE_MPI
     MPI_Init(&argc, &argv);
   #else
     static_cast<void>(argc); // elide warning about unused vars
     static_cast<void>(argv);
   #endif
   }

   ~SimpleMPIWrapper()
   {
   #ifdef TRIBOL_USE_MPI
      MPI_Finalize();
   #endif
   }
};

/*!
 * \brief Struct to hold control parameters for tests
 */
struct TestControlParameters
{
   /// Constructor
   TestControlParameters() :
        penalty_ratio            (false)
      , constant_rate_penalty    (false)
      , percent_rate_penalty     (false)
      , rate_penalty             (1.0)
      , rate_penalty_ratio       (0.0)
      , const_penalty            (1.0)
   {}

   ~TestControlParameters() 
   { 
      // no-op
   }

   real dt {0.};
   real contact_pen_frac {0.30};

   // penalty control parameters
   bool penalty_ratio;
   bool constant_rate_penalty;
   bool percent_rate_penalty;
   real rate_penalty;
   real rate_penalty_ratio;
   real const_penalty;
};

/*!
 * \brief class to hold test mesh data
 */
class TestMesh
{
public:
   /// constructor
   TestMesh();

   /// destructor
   ~TestMesh();

   /// clear function
   void clear();

   /// performs tribol registration calls and calls tribol::update()
   int tribolSetupAndUpdate( ContactMethod method,           ///< contact method
                             EnforcementMethod enforcement,  ///< constraint enforcement method
                             ContactModel model,             ///< contact model 
                             bool visualization,             ///< true if visualization
                             TestControlParameters & params  ///< control parameters struct
                            );

   /// performs tribol registration calls and calls tribol::update() using "simple" API
   int simpleTribolSetupAndUpdate( ContactMethod method,           ///< contact method
                                   EnforcementMethod enforcement,  ///< constraint enforcement method
                                   ContactModel model,             ///< contact model 
                                   bool visualization,             ///< true if visualization
                                   TestControlParameters & params  ///< control parameters struct
                            );

  /*!
   * \brief setups of a 3D contact mesh consisting of two blocks
   *
   * \param [in] numElemsX1 number of elements in the x-direction for first block
   * \param [in] numElemsY1 number of elements in the y-direction for first block
   * \param [in] numElemsZ1 number of elements in the z-direction for first block
   * \param [in] xMin1 minimum x-coordinate location of first block
   * \param [in] yMin1 minimum y-coordinate location of first block
   * \param [in] zMin1 minimum z-coordinate location of first block
   * \param [in] xMax1 maximum x-coordinate location of first block
   * \param [in] yMax1 maximum y-coordinate location of first block
   * \param [in] zMax1 maximum z-coordinate location of first block
   * \param [in] numElemsX2 number of elements in the x-direction for second block
   * \param [in] numElemsY2 number of elements in the y-direction for second block
   * \param [in] numElemsZ2 number of elements in the z-direction for second block
   * \param [in] xMin2 minimum x-coordinate location of second block
   * \param [in] yMin2 minimum y-coordinate location of second block
   * \param [in] zMin2 minimum z-coordinate location of second block
   * \param [in] xMax2 maximum x-coordinate location of second block
   * \param [in] yMax2 maximum y-coordinate location of second block
   * \param [in] zMax2 maximum z-coordinate location of second block
   * \param [in] thetaMaster angle of rotation of non-z-plane vertices about z-axis
   * \param [in] thetaSlave angle of rotation of non-z-plane vertices about z-axis
   *
   */
   void setupContactMeshHex( int numElemsX1, int numElemsY1, int numElemsZ1, 
                             real xMin1, real yMin1, real zMin1,
                             real xMax1, real yMax1, real zMax1,
                             int numElemsX2, int numElemsY2, int numElemsZ2,
                             real xMin2, real yMin2, real zMin2, 
                             real xMax2, real yMax2, real zMax2,
                             real thetaMaster, real thetaSlave );

  /*!
   * \brief sets up the Dirichlet BC node ids and values for a single 3D mesh block 
   *        for the contact PATCH TEST 
   *
   * \param [in] numElemsX number of elements in the x-direction
   * \param [in] numElemsY number of elements in the y-direction
   * \param [in] numElemsZ number of elements in the z-direction
   * \param [in] master true if this is the master block in the mesh
   * \param [in] nodeIdOffset node id offset for this block
   * \param [in] inHomogeneousGap true if enforcing gap closure through Dirichlet BCs
   * \param [in] inHomogeneousZVal z-component inhomogeneous Dirichlet BC for when inHomogeneousGap is true
   *
   * \note inHomogeneousGap is used for enforcing gap closure with Dirichlet BCs  instead of
   *       contact enforcement. This is used in tribol/tests/tribol_mortar_pressure_sol.cpp
   *
   */
   void setupPatchTestDirichletBCs( int numElemsX, int numElemsY, int numElemsZ, 
                                    bool master, int nodeIdOffset, 
                                    bool inHomogeneousGap, 
                                    real inHomogeneousZVal = 0. );

  /*!
   * \brief sets up pressure dof ids for a 3D slave mesh block for PATCH TEST
   *
   * \param [in] numElemsX number of elements in x-direction of slave block
   * \param [in] numElemsY number of elements in y-direction of slave block
   * \param [in] numElemsZ number of elements in z-direction of slave block
   * \param [in] nodeIdOffset slave node id offset
   * \param [in] contact true if enforcing zero gap using contact enforcement
   * \param [in,out] presDofs pointer to array of slave pressure dof node ids
   *
   */
   void setupPatchTestPressureDofs( int numElemsX, int numElemsY, int numElemsZ, 
                                    int nodeIdOffset, bool contact, bool master );

  /*!
   * \brief sets up an mfem mesh object representation of the original test mesh
   *
   */
   void setupMfemMesh( );

  /*!
   * \brief allocates and sets velocity arrays
   *  
   * \param [in] valX x-velocity value to set velocity arrays
   * \param [in] valY y-velocity value to set velocity arrays
   * \param [in] valZ z-velocity value to set velocity arrays
   *
   * \note in the future we should handle an array of velocity values
   *
   */
   void allocateAndSetVelocities( int meshId, real valX, real valY, real valZ=0. );

  /*!
   * \brief allocates and sets bulk modulus arrays on mesh
   *  
   * \param [in] val single scalar bulk modulus 
   *
   * \note in the future we should handle an array of values
   *
   */
   void allocateAndSetBulkModulus( int meshId, real val );

  /*!
   * \brief allocates and sets element thickness arrays on mesh
   *  
   * \param [in] val single scalar element thickness
   *
   * \note this is suitable for a uniform mesh. In the future we 
   *       should handle an array of values
   *
   */
   void allocateAndSetElementThickness( int meshId, real t );

  /*!
   * \brief wraps element Jacobian calculations for linear elasticity
   *
   * \param [in,out] A      mfem sparse matrix data object where Jac. contributions are summed into
   * \param [in]     nu     Poisson's ratio for linear elastic material
   * \param [in]     youngs Young's modulus for linear elastic material
   *
   *
   * \pre TestMesh::mfem_mesh != nullptr
   * \pre matrix A cannot have been finalized yet
   */
   void computeEquilibriumJacobian( mfem::SparseMatrix * const A, 
                                    real const nu, real const youngs );

  /*!
   * \brief wraps element Jacobian calculations for linear elasticity 
   *
   * \param [in,out] A      mfem sparse matrix data object where Jac. contributions are summed into
   * \param [in]     eInteg mfem elasticity integrator
   * \param [in]     fes    finite element space
   *
   *
   * \pre fes != nullptr
   * \pre TestMesh::mfem_mesh != nullptr
   * \pre matrix A cannot have been finalized yet
   */
   void computeEquilibriumJacobian( mfem::SparseMatrix * const A, 
                                    mfem::ElasticityIntegrator * eInteg,
                                    mfem::FiniteElementSpace * fes );

  /*!
   * \brief computes all element Jacobian contributions for linear elasticity using MFEM
   *
   * \param [in,out] A      mfem sparse matrix data object where Jac. contributions are summed into
   * \param [in]     eInt   mfem elasticity integrator
   * \param [in]     fes    finite element space
   *
   *
   * \pre fes != nullptr
   * \pre TestMesh::mfem_mesh != nullptr
   * \pre matrix A cannot have been finalized yet
   */
   void computeElementJacobianContributions( mfem::SparseMatrix * const A,
                                             mfem::ElasticityIntegrator * eInt,
                                             mfem::FiniteElementSpace * fe_space,
                                             bool matrixDebug = false );

  /*!
   * \brief takes an oversized tribol sparse matrix and condenses it to an appropriate
   *        system size sparse matrix
   *
   * \param [in,out] ATribol Tribol-sized dense matrix (dim * totalNumberOfNodes + totalNumberOfNodes)
   * \param [in,out] ASystem system-sized sparse matrix (dim * totalNumberOfNodes + numSlaveNodes)
   *
   * \pre ATribol != nullptr
   * \pre ASystem != nullptr
   *
   * \note this routine is a little odd taking a dense matrix and populating a sparse
   *       matrix. This reflects the way the Lagrange multiplier patch test code 
   *       was originally written and does not represent the optimal way of doing things
   */
   void tribolMatrixToSystemMatrix( mfem::DenseMatrix * const ATribol,
                                    mfem::SparseMatrix * const ASystem );

  /*!
   * \brief populate a right hand side length vector with gap evaluations
   *
   * \param [in,out] v pointer to right hand side vector
   *
   * \pre v of length, dim * numTotalNodes + numPressureDofs
   *
   */
   void getGapEvals( real * const v );

  /*!
   * \brief Modifies matrix A and rhs vector b to enforce Dirichlet BCs
   *
   * \param [in,out] A mfem sparse matrix
   * \param [in,out] b mfem vector (rhs vector)
   * \param [in]     contact applies BCs different for a contact vs. non-contact solution
   *
   * \note The boolean, contact, is something utilized for the SINGLE_MORTAR LM 
   *       patch unit test
   *
   * \pre b of length, dim * numTotalNodes + numPressureDofs
   * \pre A of dimensions (length b) x (length b)
   *
   */
   void enforceDirichletBCs( mfem::SparseMatrix * const A,
                             mfem::Vector * const b,
                             bool contact = true );

   /// print mesh to vtk file
   void testMeshToVtk( const std::string& dir, ///< Name of the output directory
                       int cycle,              ///< Cycle number
                       double time             ///< Simulation time 
                     );

public:

   mfem::Mesh * mfem_mesh; 

   // public member variables
   int masterMeshId;         ///< Mesh id for master portion of mesh
   int slaveMeshId;          ///< Mesh id for slave portion of mesh
   int numTotalNodes;        ///< Total number of nodes in the mesh 
   int numMasterNodes;       ///< Number of master nodes (not just surface nodes)
   int numSlaveNodes;        ///< Number of slave nodes (not just surface nodes)
   int numSlaveSurfaceNodes; ///< Number of surface nodes on the slave side
   int numTotalElements;     ///< Total number of elements 
   int numMasterElements;    ///< Number of master elements
   int numSlaveElements;     ///< Number of slave nodes 
   int numTotalFaces;        ///< Total number of faces
   int numMasterFaces;       ///< Number of master faces
   int numSlaveFaces;        ///< Number of slave faces 
   int numNodesPerFace;      ///< Number of nodes per face
   int numNodesPerElement;   ///< Number of nodes per element
   int dim;                  ///< Mesh dimension
   int *dirNodesX1;          ///< Pointer to list of master node ids with x-component Dirichlet BCs 
   int *dirNodesY1;          ///< Pointer to list of master node ids with y-component Dirichlet BCs
   int *dirNodesZ1;          ///< Pointer to list of master node ids with z-component Dirichlet BCs 
   double *iDirValX1;        ///< Pointer to x-component Dirichlet BC values for specified master nodes
   double *iDirValY1;        ///< Pointer to y-component Dirichlet BC values for specified master nodes
   double *iDirValZ1;        ///< Pointer to z-component Dirichlet BC values for specified master nodes
   int *dirNodesX2;          ///< Pointer to list of slave node ids with x-component Dirichlet BCs
   int *dirNodesY2;          ///< Pointer to list of slave node ids with y-component Dirichlet BCs
   int *dirNodesZ2;          ///< Pointer to list of slave node ids with z-component Dirichlet BCs
   double *iDirValX2;        ///< Pointer to x-component Dirichlet BC values for specified slave nodes
   double *iDirValY2;        ///< Pointer to y-component Dirichlet BC values for specified slave nodes
   double *iDirValZ2;        ///< Pointer to z-component Dirichlet BC values for specified slave nodes
   int *faceConn1;           ///< Pointer to master face connectivity
   int *faceConn2;           ///< Pointer to slave face connectivity
   int *elConn1;             ///< Pointer to master element connectivity
   int *elConn2;             ///< Pointer to slave element connectivity 
   int *presDofs1;           ///< Pointer to master node ids with a pressure BC
   int *presDofs2;           ///< Pointer to slave node ids with a pressure BC

   double *fx1; ///< Master nodal forces, x-component 
   double *fy1; ///< Master nodal forces, y-component 
   double *fz1; ///< Master nodal forces, z-component
   double *fx2; ///< Slave nodal forces, x-component
   double *fy2; ///< Slave nodal forces, y-component
   double *fz2; ///< Slave nodal forces, z-component
   double *vx1; ///< Master nodal velocities, x-component
   double *vy1; ///< Master nodal velocities, y-component
   double *vz1; ///< Master nodal velocities, z-component
   double *vx2; ///< Slave nodal velocities, x-component
   double *vy2; ///< Slave nodal velocities, y-component
   double *vz2; ///< Slave nodal velocities, z-component
   double *x;   ///< Nodal coordinates, x-component
   double *y;   ///< Nodal coordinates, y-component
   double *z;   ///< Nodal coordinates, z-component

   // CSR storage
   int* I;       ///< Offsets for CSR storage
   int* J;       ///< Column indices for all nonzero entries in CSR storage
   double *vals; ///< Array of nonzero values 

   double *gaps;      ///< Array of nodal gaps
   double *pressures; ///< Array of nodal pressures

   real * master_bulk_mod;
   real * master_element_thickness;
   real * slave_bulk_mod;
   real * slave_element_thickness;

   bool registered_velocities1;
   bool registered_velocities2;

  public:
   double* getX() const {return x;}
   double* getY() const {return y;}
   double* getZ() const {return z;}
   int* getMasterFaceConnectivity() const {return faceConn1;}
   int* getSlaveFaceConnectivity()  const {return faceConn2;}

   int getNumTotalNodes() const { return numTotalNodes;}

  public:
   /// Needed for initial shroud interface
   int getMasterFaceConnectivitySize() const { return numNodesPerFace * numMasterFaces ; }
   int getSlaveFaceConnectivitySize() const { return numNodesPerFace * numSlaveFaces ; }

protected:
   // NONE
};

} // end of namespace "tribol"

#endif /* SRC_UTILS_TESTUTILS_HPP_ */
