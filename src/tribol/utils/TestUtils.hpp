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
   * \param [in] thetaMortar angle of rotation of non-z-plane vertices about z-axis
   * \param [in] thetaNonmortar angle of rotation of non-z-plane vertices about z-axis
   *
   */
   void setupContactMeshHex( int numElemsX1, int numElemsY1, int numElemsZ1, 
                             real xMin1, real yMin1, real zMin1,
                             real xMax1, real yMax1, real zMax1,
                             int numElemsX2, int numElemsY2, int numElemsZ2,
                             real xMin2, real yMin2, real zMin2, 
                             real xMax2, real yMax2, real zMax2,
                             real thetaMortar, real thetaNonmortar );

  /*!
   * \brief sets up the Dirichlet BC node ids and values for a single 3D mesh block 
   *        for the contact PATCH TEST 
   *
   * \param [in] numElemsX number of elements in the x-direction
   * \param [in] numElemsY number of elements in the y-direction
   * \param [in] numElemsZ number of elements in the z-direction
   * \param [in] mortar true if this is the mortar block in the mesh
   * \param [in] nodeIdOffset node id offset for this block
   * \param [in] inHomogeneousGap true if enforcing gap closure through Dirichlet BCs
   * \param [in] inHomogeneousZVal z-component inhomogeneous Dirichlet BC for when inHomogeneousGap is true
   *
   * \note inHomogeneousGap is used for enforcing gap closure with Dirichlet BCs  instead of
   *       contact enforcement. This is used in tribol/tests/tribol_mortar_pressure_sol.cpp
   *
   */
   void setupPatchTestDirichletBCs( int numElemsX, int numElemsY, int numElemsZ, 
                                    bool mortar, int nodeIdOffset, 
                                    bool inHomogeneousGap, 
                                    real inHomogeneousZVal = 0. );

  /*!
   * \brief sets up pressure dof ids for a 3D nonmortar mesh block for PATCH TEST
   *
   * \param [in] numElemsX number of elements in x-direction of nonmortar block
   * \param [in] numElemsY number of elements in y-direction of nonmortar block
   * \param [in] numElemsZ number of elements in z-direction of nonmortar block
   * \param [in] nodeIdOffset nonmortar node id offset
   * \param [in] contact true if enforcing zero gap using contact enforcement
   * \param [in,out] presDofs pointer to array of nonmortar pressure dof node ids
   *
   */
   void setupPatchTestPressureDofs( int numElemsX, int numElemsY, int numElemsZ, 
                                    int nodeIdOffset, bool contact, bool mortar );

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
   * \param [in,out] ASystem system-sized sparse matrix (dim * totalNumberOfNodes + numNonmortarNodes)
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
   int mortarMeshId;         ///< Mesh id for mortar portion of mesh
   int nonmortarMeshId;          ///< Mesh id for nonmortar portion of mesh
   int numTotalNodes;        ///< Total number of nodes in the mesh 
   int numMortarNodes;       ///< Number of mortar nodes (not just surface nodes)
   int numNonmortarNodes;        ///< Number of nonmortar nodes (not just surface nodes)
   int numNonmortarSurfaceNodes; ///< Number of surface nodes on the nonmortar side
   int numTotalElements;     ///< Total number of elements 
   int numMortarElements;    ///< Number of mortar elements
   int numNonmortarElements;     ///< Number of nonmortar nodes 
   int numTotalFaces;        ///< Total number of faces
   int numMortarFaces;       ///< Number of mortar faces
   int numNonmortarFaces;        ///< Number of nonmortar faces 
   int numNodesPerFace;      ///< Number of nodes per face
   int numNodesPerElement;   ///< Number of nodes per element
   int dim;                  ///< Mesh dimension
   int *dirNodesX1;          ///< Pointer to list of mortar node ids with x-component Dirichlet BCs 
   int *dirNodesY1;          ///< Pointer to list of mortar node ids with y-component Dirichlet BCs
   int *dirNodesZ1;          ///< Pointer to list of mortar node ids with z-component Dirichlet BCs 
   double *iDirValX1;        ///< Pointer to x-component Dirichlet BC values for specified mortar nodes
   double *iDirValY1;        ///< Pointer to y-component Dirichlet BC values for specified mortar nodes
   double *iDirValZ1;        ///< Pointer to z-component Dirichlet BC values for specified mortar nodes
   int *dirNodesX2;          ///< Pointer to list of nonmortar node ids with x-component Dirichlet BCs
   int *dirNodesY2;          ///< Pointer to list of nonmortar node ids with y-component Dirichlet BCs
   int *dirNodesZ2;          ///< Pointer to list of nonmortar node ids with z-component Dirichlet BCs
   double *iDirValX2;        ///< Pointer to x-component Dirichlet BC values for specified nonmortar nodes
   double *iDirValY2;        ///< Pointer to y-component Dirichlet BC values for specified nonmortar nodes
   double *iDirValZ2;        ///< Pointer to z-component Dirichlet BC values for specified nonmortar nodes
   int *faceConn1;           ///< Pointer to mortar face connectivity
   int *faceConn2;           ///< Pointer to nonmortar face connectivity
   int *elConn1;             ///< Pointer to mortar element connectivity
   int *elConn2;             ///< Pointer to nonmortar element connectivity 
   int *presDofs1;           ///< Pointer to mortar node ids with a pressure BC
   int *presDofs2;           ///< Pointer to nonmortar node ids with a pressure BC

   double *fx1; ///< Mortar nodal forces, x-component 
   double *fy1; ///< Mortar nodal forces, y-component 
   double *fz1; ///< Mortar nodal forces, z-component
   double *fx2; ///< Nonmortar nodal forces, x-component
   double *fy2; ///< Nonmortar nodal forces, y-component
   double *fz2; ///< Nonmortar nodal forces, z-component
   double *vx1; ///< Mortar nodal velocities, x-component
   double *vy1; ///< Mortar nodal velocities, y-component
   double *vz1; ///< Mortar nodal velocities, z-component
   double *vx2; ///< Nonmortar nodal velocities, x-component
   double *vy2; ///< Nonmortar nodal velocities, y-component
   double *vz2; ///< Nonmortar nodal velocities, z-component
   double *x;   ///< Nodal coordinates, x-component
   double *y;   ///< Nodal coordinates, y-component
   double *z;   ///< Nodal coordinates, z-component

   // CSR storage
   int* I;       ///< Offsets for CSR storage
   int* J;       ///< Column indices for all nonzero entries in CSR storage
   double *vals; ///< Array of nonzero values 

   double *gaps;      ///< Array of nodal gaps
   double *pressures; ///< Array of nodal pressures

   real * mortar_bulk_mod;
   real * mortar_element_thickness;
   real * nonmortar_bulk_mod;
   real * nonmortar_element_thickness;

   bool registered_velocities1;
   bool registered_velocities2;

  public:
   double* getX() const {return x;}
   double* getY() const {return y;}
   double* getZ() const {return z;}
   int* getMortarFaceConnectivity() const {return faceConn1;}
   int* getNonmortarFaceConnectivity()  const {return faceConn2;}

   int getNumTotalNodes() const { return numTotalNodes;}

  public:
   /// Needed for initial shroud interface
   int getMortarFaceConnectivitySize() const { return numNodesPerFace * numMortarFaces ; }
   int getNonmortarFaceConnectivitySize() const { return numNodesPerFace * numNonmortarFaces ; }

protected:
   // NONE
};

} // end of namespace "tribol"

namespace mfem_ext
{

/// Simple central difference method with homogeneous velocity BCs.
class CentralDiffSolver : public mfem::SecondOrderODESolver
{
public:
   /**
    * @brief Construct a new central difference solver object
    *
    * Supports homogeneous velocity BCs
    * 
    * @param bc_vdofs_ List of vdofs to set homogeneous velocity BCs
    */
   CentralDiffSolver(const mfem::Array<int>& bc_vdofs_);

   /**
    * @brief Updates x and dxdt after taking a step of size dt
    * 
    * @param x Displacement vector
    * @param dxdt Velocity vector
    * @param t Current time
    * @param dt Timestep size
    */
   void Step(mfem::Vector& x, mfem::Vector& dxdt, double& t, double& dt) override;

private:
   /**
    * @brief Acceleration vector
    */
   mfem::Vector accel;

   /**
    * @brief List of vdofs to apply homogeneous velocity BCs
    */
   mfem::Array<int> bc_vdofs;

   /**
    * @brief Tracks whether a step has been taken yet
    */
   bool first_step;

   /**
    * @brief Applies homogeneous BCs to dxdt
    * 
    * @param dxdt Velocity vector
    */
   void SetHomogeneousBC(mfem::Vector& dxdt) const;
};

/// Explicit solid mechanics update (lumped mass, optional damping)
class ExplicitMechanics : public mfem::SecondOrderTimeDependentOperator
{
public:
   /**
    * @brief Construct a new explicit mechanics operator (elasticity, lumped
    * mass, optional damping, and external force vector)
    *
    * @param fespace Finite element space of the primary fields
    * @param rho Density
    * @param lambda Lame constant
    * @param mu Lame constant
    * @param damping_ Damping matrix
    */
   ExplicitMechanics(
      mfem::ParFiniteElementSpace& fespace, 
      mfem::Coefficient& rho,
      mfem::Coefficient& lambda,
      mfem::Coefficient& mu
   );

   /**
    * @brief Compute acceleration given displacement and velocity
    * 
    * @param u Displacement vector
    * @param dudt Velocity vector
    * @param a Acceleration vector
    */
   void Mult(
      const mfem::Vector& u,
      const mfem::Vector& dudt,
      mfem::Vector& a
   ) const override;

   /**
    * @brief External force contribution (must be manually updated)
    */
   mfem::Vector f_ext;

private:
   /**
    * @brief Elasticity bilinear form
    */
   mfem::ParBilinearForm elasticity;
   
   /**
    * @brief Inverse lumped mass matrix
    */
   mfem::Vector inv_lumped_mass;
};

} // end of namespace "mfem_ext"

#endif /* SRC_UTILS_TESTUTILS_HPP_ */
