// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Tribol Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (MIT)

// Tribol includes
#include "tribol/types.hpp"
#include "tribol/interface/tribol.hpp"
#include "tribol/utils/TestUtils.hpp"
#include "tribol/common/Parameters.hpp"
#include "tribol/mesh/MethodCouplingData.hpp"
#include "tribol/mesh/CouplingSchemeManager.hpp"
#include "tribol/mesh/CouplingScheme.hpp"
#include "tribol/mesh/MeshData.hpp"
#include "tribol/mesh/MeshManager.hpp"
#include "tribol/physics/Mortar.hpp"
#include "tribol/physics/AlignedMortar.hpp"
#include "tribol/geom/GeomUtilities.hpp"

// Axom includes
#include "axom/slic.hpp"

// MFEM includes
#include "mfem.hpp"

// gtest includes
#include "gtest/gtest.h"

// c++ includes
#include <cmath> // std::abs, std::cos, std::sin
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>

using real = tribol::real;

/*!
 * Test fixture class with some setup necessary to test
 * the Lagrange multiplier pressure field solution for 
 * Tribol's simplified mortar method on a single element 
 * on single element contact problem
 */
class MortarLMPatchTest : public ::testing::Test
{
   
public:
   tribol::TestMesh m_mesh;

   void computeContactSolution( int nMortarElemsX, int nMortarElemsY, int nMortarElemsZ,
                                int nNonmortarElemsX, int nNonmortarElemsY, int nNonmortarElemsZ,
                                real xMin1, real yMin1, real zMin1,
                                real xMax1, real yMax1, real zMax1,
                                real xMin2, real yMin2, real zMin2,
                                real xMax2, real yMax2, real zMax2,
                                real thetaM, real thetaS,
                                tribol::ContactMethod method,
                                std::string test_name,
                                bool cntct, bool writeOutput, bool debug, bool vis,
                                mfem::Vector &sol, real &pressure_rel_error );

protected:

   void SetUp() override
   {
   }

   void TearDown() override
   {
      // call clear() on mesh object to be safe
      this->m_mesh.clear();
   }

protected:

};

void MortarLMPatchTest::computeContactSolution( int nMortarElemsX, int nMortarElemsY, int nMortarElemsZ,
                                                int nNonmortarElemsX, int nNonmortarElemsY, int nNonmortarElemsZ,
                                                real xMin1, real yMin1, real zMin1,
                                                real xMax1, real yMax1, real zMax1,
                                                real xMin2, real yMin2, real zMin2,
                                                real xMax2, real yMax2, real zMax2,
                                                real thetaM, real thetaS,
                                                tribol::ContactMethod method,
                                                std::string test_name,
                                                bool cntct, bool writeOutput, bool TRIBOL_UNUSED_PARAM(debug), bool vis,
                                                mfem::Vector &sol, real &pressure_rel_error )
{
   bool inHomogeneous = false;
   real inHomogeneousVal = 0.;
   if (!cntct)
   {
      inHomogeneous = true;
      inHomogeneousVal = 0.05; // hard coded for this example
   }

   this->m_mesh.setupContactMeshHex( nMortarElemsX, nMortarElemsY, nMortarElemsZ, 
                                     xMin1, yMin1, zMin1, xMax1, yMax1, zMax1,
                                     nNonmortarElemsX, nNonmortarElemsY, nNonmortarElemsZ, 
                                     xMin2, yMin2, zMin2, xMax2, yMax2, zMax2,
                                     thetaM, thetaS );

   // setup mortar boundary conditions
   this->m_mesh.setupPatchTestDirichletBCs( this->m_mesh.mortarMeshId, nMortarElemsX, nMortarElemsY, nMortarElemsZ, 
                                            0, inHomogeneous, -inHomogeneousVal );

   // setup nonmortar boundary conditions
   this->m_mesh.setupPatchTestDirichletBCs( this->m_mesh.nonmortarMeshId, nNonmortarElemsX, nNonmortarElemsY, nNonmortarElemsZ, 
                                            this->m_mesh.numMortarNodes, inHomogeneous, inHomogeneousVal );

   // setup DUMMY MORTAR pressure dof array. Consider getting rid of 
   // mortar pressure dofs based on new way stiffness data is handled 
   // after passed back from Tribol. ALWAYS call this function with 
   // these arguments for the mortar block
   this->m_mesh.setupPatchTestPressureDofs( this->m_mesh.mortarMeshId, nMortarElemsX, nMortarElemsY, nMortarElemsZ,
                                            0, false );

   // setup NONMORTAR pressure dofs
   this->m_mesh.setupPatchTestPressureDofs( this->m_mesh.nonmortarMeshId, nNonmortarElemsX, nNonmortarElemsY, nNonmortarElemsZ,
                                            this->m_mesh.numMortarNodes, true );

   // specify if contact is on
   //bool matrixDebug = debug; // this controls whether the matrix printed is without (true) or with BCs applied 
   bool output = writeOutput;
   bool contact = cntct;

   // register mesh, fields and coupling scheme with Tribol and call update. Note, the mfem 
   // mesh is not needed here. The mfem mesh is simply for the equilibrium calculations
   tribol::TestControlParameters parameters;

   int test_mesh_update_err = 
      this->m_mesh.tribolSetupAndUpdate( method, tribol::LAGRANGE_MULTIPLIER, 
                                         tribol::FRICTIONLESS, vis, parameters );

   EXPECT_EQ( test_mesh_update_err, 0 );

   // setup mfem mesh
   this->m_mesh.setupMfemMesh( );

   // setup a temporary stacked array of nodal coordinates for use in the instantiation 
   // of a reference configuration grid function. The nodal coordinates are stacked x, then 
   // y, then z. Also setup a local interleaved array of incremental nodal displacements.
   real xyz[ this->m_mesh.dim * this->m_mesh.numTotalNodes ];
   real xyz_inc[ this->m_mesh.dim * this->m_mesh.numTotalNodes ];
   for (int i=0; i<this->m_mesh.numTotalNodes; ++i)
   {
      xyz[ i ] = this->m_mesh.x[i];
      xyz[ this->m_mesh.numTotalNodes + i ]   = this->m_mesh.y[i];
      xyz[ 2*this->m_mesh.numTotalNodes + i ] = this->m_mesh.z[i];

      xyz_inc[ this->m_mesh.dim * i ]     = 0.;
      xyz_inc[ this->m_mesh.dim * i + 1 ] = 0.;
      xyz_inc[ this->m_mesh.dim * i + 2 ] = 0.;
   }

   // define the FE collection and finite element space
   mfem::FiniteElementSpace * fe_space { nullptr };
   int order = 1;
   mfem::H1_FECollection fe_coll( order, this->m_mesh.dim );
   fe_space = new mfem::FiniteElementSpace( this->m_mesh.mfem_mesh, &fe_coll, this->m_mesh.dim );

   // get Jacobian sparse matrix from Tribol
   mfem::SparseMatrix * tribolJac { nullptr };
   int sparseMatErr = tribol::getMfemSparseMatrix( &tribolJac, 0 );

   EXPECT_EQ( sparseMatErr, 0 ); 

   // add hex8 equilibrium contributions for linear elasticity
   // define the lambda and mu constant coefficients
   double lambda_val, mu_val, nu_val, youngs_val;
   nu_val = 0.33;
   youngs_val = 3.E7; 
   lambda_val = (youngs_val * nu_val) / ((1. + nu_val) * (1. - 2. * nu_val));
   mu_val = youngs_val / (2.*(1.+nu_val));
   mfem::ConstantCoefficient lambda(lambda_val);
   mfem::ConstantCoefficient mu(mu_val);

   // instantiate elasticity integrator
   mfem::ElasticityIntegrator elastInteg( lambda, mu );

   // compute equilibrium contributions and sum into Jacobian
   m_mesh.computeEquilibriumJacobian( tribolJac, &elastInteg, fe_space );

   // make sure sparse matrix is finalized to convert to CSR format
   tribolJac->Finalize();

   // convert sparse matrix to dense matrix for output
   mfem::DenseMatrix dTribolJac;
   tribolJac->ToDenseMatrix(dTribolJac);

   // instantiate mfem vector for rhs vector. The length of this vector is the 
   // (space dimension) x (total number of mesh nodes) + (number of nonmortar nodes in contact), 
   // where the last addition is for the pressure lagrange multiplier field
   int rhs_size = this->m_mesh.dim * this->m_mesh.numTotalNodes +
                  this->m_mesh.numNonmortarSurfaceNodes;
   double b[ rhs_size ];

   // initialize b vector
   for (int i=0; i<rhs_size; ++i)
   {
      b[i] = 0.;
   }

   // instantiate mfem vector for right hand side
   mfem::Vector rhs( &b[0], rhs_size );

   ////////////////////////////////////////////////////
   //                                                //
   // Condense the Tribol Jacobian into a system     //
   // that is (dim*numTotalNodes + numPressureDof).  //
   // Place contributions into another sparse matrix //
   // sized accordingly.                             //
   //                                                //
   ////////////////////////////////////////////////////

   int solveSize = rhs_size;
   mfem::SparseMatrix jac( solveSize );

   int numRows = jac.NumRows();

   this->m_mesh.tribolMatrixToSystemMatrix( &dTribolJac, &jac );

   // note: this does not populate the right hand side with any contact weak form 
   // residual terms. That is, the initial guess for pressure is zero, and therefore
   // there are no contact contributions to the equilibrium residual. Here we only 
   // assemble gap contributions. Recall at this point that tribol::update() is 
   // called again with the updated contact solution, so doing so would actually 
   // mess with the RHS in a negative way
   if (contact)
   {
      this->m_mesh.getGapEvals( &b[0] );
   }

   // zero out all Dirichlet BC components for each block
   this->m_mesh.enforceDirichletBCs( &jac, &rhs, contact );

   jac.Finalize();
   mfem::DenseMatrix dJac;
   jac.ToDenseMatrix(dJac);
   int rank = dJac.Rank(1.e-15);

   SLIC_INFO( "Matrix rank: " << rank );

   if (rank < numRows)
   {
      SLIC_ERROR("Jacobian rank (" << rank << ") less than row dimension (" << numRows << ")" );
   }

   // instantiate mfem dense matrix inverse object and 
   // solution vector
   mfem::DenseMatrixInverse invJ( dJac );

   // Solve the system
   invJ.Mult( rhs, sol );

   // get solution data
   double * sol_data = sol.GetData();

   // update the tribol nodal pressure array. This traverses 
   // the connectivity array, which will populate nonmortar pressures 
   // more than once for nonmortar nodes shared between multiple faces, 
   // but this is easy for now (SRW). To do this we need a mapping 
   // from local nonmortar pressure dofs (in the solution vector) to 
   // the global mesh connectivity ids. I don't create this mapping, 
   // instead I exploit the fact that the nonmortar pressure nodes are 
   // ordered consecutively with an offset that is the number of 
   // mortar nodes in the mortar block
   double nonmortarForceSum = 0.;
   if (contact)
   {
      int connSize = this->m_mesh.numNonmortarFaces * this->m_mesh.numNodesPerFace;
      for (int i=0; i<connSize; ++i)
      {
         int offset = this->m_mesh.dim * (this->m_mesh.numMortarNodes + this->m_mesh.numNonmortarNodes);
         int nonmortarOffset = this->m_mesh.numMortarNodes;
         int id = this->m_mesh.faceConn2[i];
         this->m_mesh.pressures[id] = sol_data[ offset + id - nonmortarOffset ];
      }

      // zero out mortar and nonmortar nodal force contributions for equilibrium 
      // residual evaluation. Note, don't update the nodal coordinates! We want 
      // to evaluate the equilibrium residual with the current pressure 
      // solution and the contact overlaps (i.e. mortar weights) used to 
      // solve for those pressures. The updated nodal coordinates should be 
      // used for a gap only evaluation, which is currently not done.
      for (int i=0; i<this->m_mesh.numTotalNodes; ++i)
      {
         this->m_mesh.fx1[i] = 0.;
         this->m_mesh.fy1[i] = 0.;
         this->m_mesh.fz1[i] = 0.;
         this->m_mesh.fx2[i] = 0.;
         this->m_mesh.fy2[i] = 0.;
         this->m_mesh.fz2[i] = 0.;
      }

      // call tribol update() again for residual only evaluation
      // TODO check to make sure this call here after refactor works, SRW
      tribol::setLagrangeMultiplierOptions( 0, tribol::ImplicitEvalMode::MORTAR_RESIDUAL, 
                                            tribol::SparseMode::MFEM_LINKED_LIST );

      double dt = 1.0;
      int tribol_update_err = tribol::update( 1, 1., dt );

      EXPECT_EQ( tribol_update_err, 0 );

      // sum the nodal force contributions on the nonmortar side
      for (int i=0; i<this->m_mesh.numNonmortarSurfaceNodes; ++i)
      {
         int offset = this->m_mesh.numMortarNodes;
         // exploit offset and contiguous node numbering in the indexing here.
         nonmortarForceSum += this->m_mesh.fz2[ offset + i ];
      }
      pressure_rel_error = nonmortarForceSum;
      SLIC_INFO("NODAL FORCE SUM (NONMORTAR, TRIBOL RESIDUALS): " << nonmortarForceSum);
   }

   // update nodal coordinates in separate stacked array. Keep original 
   // mesh array as that of the reference configuration
   for (int i=0; i<this->m_mesh.numTotalNodes; ++i)
   {
      for (int j=0; j<this->m_mesh.dim; ++j)
      { 
         xyz[ this->m_mesh.numTotalNodes * j + i ] 
            += sol_data[ this->m_mesh.dim * i + j ];
         xyz_inc[ this->m_mesh.numTotalNodes * j + i ] 
            = sol_data[ this->m_mesh.dim * i + j ];
      }
   }

   // compute stress update
   mfem::GridFunction u( fe_space, &xyz_inc[0] );
   const int tdim = this->m_mesh.dim*(this->m_mesh.dim+1)/2;
   mfem::FiniteElementSpace flux_fespace( this->m_mesh.mfem_mesh, &fe_coll, tdim );
   mfem::GridFunction stress( &flux_fespace );
   stress = 0.;
   mfem::BilinearFormIntegrator *integ = &elastInteg;

   u.ComputeFlux( *integ, stress );

   // instantiate stress33 relative error grid function
   mfem::GridFunction stress33_rel_error( stress );

   // compute pressure difference from computed contact pressure solution
   // and stress solution from mfem
   if (contact)
   {
      mfem::Vector node_vals;
      stress.GetNodalValues( node_vals, 3 );
      real *data = stress.GetData();
      pressure_rel_error 
         -= data[2*this->m_mesh.numTotalNodes + this->m_mesh.faceConn2[0]];
      pressure_rel_error 
         /= data[2*this->m_mesh.numTotalNodes + this->m_mesh.faceConn2[0]];
      pressure_rel_error = std::abs(pressure_rel_error);

      real * error_data = stress33_rel_error.GetData();
      for (int i=(2*this->m_mesh.numTotalNodes); i<(2*this->m_mesh.numTotalNodes+this->m_mesh.numTotalNodes); ++i)
      {
         error_data[i] -= nonmortarForceSum;
         error_data[i] /= nonmortarForceSum;
         error_data[i] = std::abs(error_data[i]);
      }
   }

   if (vis)
   {
      // mesh output
      std::ostringstream mesh_ref_name, mesh_cur_name;
      mesh_ref_name << "mesh_ref_" << test_name << "." << std::setfill('0') << std::setw(6) << ".vtk";
      mesh_cur_name << "mesh_cur_" << test_name << "." << std::setfill('0') << std::setw(6) << ".vtk";

      std::ofstream mesh_ref_ofs(mesh_ref_name.str().c_str());
      mesh_ref_ofs.precision(8);
 
      // print reference configuration mesh
      this->m_mesh.mfem_mesh->PrintVTK( mesh_ref_ofs, 0, 0 );

      // set the current configuration vector and mesh node grid function for output
      mfem::Vector x_cur( &xyz[0], this->m_mesh.dim*this->m_mesh.numTotalNodes );
      this->m_mesh.mfem_mesh->SetNodes( x_cur );

      // print current mesh
      std::ofstream mesh_cur_ofs(mesh_cur_name.str().c_str());
      mesh_cur_ofs.precision(8);
      this->m_mesh.mfem_mesh->PrintVTK( mesh_cur_ofs, 0, 0 );

      // save stress in current mesh .vtk. Note, do this for contact and no contact
      stress.SaveVTK( mesh_cur_ofs, "stress", 0 );
      stress33_rel_error.SaveVTK( mesh_cur_ofs, "stress_rel_error", 0 );

   }

   // DEBUG: print matrix and vector output. Leave this code in here 
   // and guard with boolean
   if (output)
   {
      std::ofstream matrix;
      matrix.setf(std::ios::scientific);
      matrix.precision(0);
      std::ostringstream suffix_matrix;
      suffix_matrix << "jacobian_" << test_name << ".txt";
      matrix.open("matrix_" + suffix_matrix.str());

      std::ofstream sol_vec;
      sol_vec.setf(std::ios::scientific);
      sol_vec.precision(2);
      std::ostringstream suffix_sol;
      suffix_sol << "vector_" << test_name << ".txt";
      sol_vec.open("sol_" + suffix_sol.str());

      std::ofstream rhs_vec;
      rhs_vec.setf(std::ios::scientific);
      rhs_vec.precision(2);
      std::ostringstream suffix_rhs;
      suffix_rhs << "vector_" << test_name << ".txt";
      rhs_vec.open("rhs_" + suffix_rhs.str());

      for (int i=0; i<jac.NumRows(); ++i)
      {
         rhs_vec << b[i] << "\n";
         sol_vec << sol_data[i] << "\n";
         for (int j=0; j<jac.NumCols(); ++j)
         {
            double val = dJac(i,j);
            matrix << val << "  ";
         }
         matrix << "\n";
      }

      sol_vec.close();
      rhs_vec.close();
      matrix.close();
   }

   tribol::finalize();

   // delete the jacobian matrix
   delete fe_space;

} // end computeContactSolution

TEST_F( MortarLMPatchTest, single_mortar_uniform_patch )
{

   mfem::Vector xsol_base;
   mfem::Vector xsol_cntct;

   int nMortarElems = 4;
   int nElemsXM = nMortarElems;
   int nElemsYM = nMortarElems;
   int nElemsZM = 1; 

   int nNonmortarElems = 4;
   int nElemsXS = nNonmortarElems;
   int nElemsYS = nNonmortarElems;
   int nElemsZS = 1;

   real x_min1 = 0.;
   real y_min1 = 0.;
   real z_min1 = 0.; 
   real x_max1 = 1.;
   real y_max1 = 1.;
   real z_max1 = 1.05;

   real x_min2 = 0.;
   real y_min2 = 0.;
   real z_min2 = 0.95;
   real x_max2 = 1.;
   real y_max2 = 1.;
   real z_max2 = 2.;

   bool output = false; // solution vector and rhs
   bool contact = true;
   bool debug = false;  // element matrix writes to .txt files
   bool visualization = false; // output of .vtk with mesh and stress output

   real thetaM = 0.;
   real thetaS = 0.;

   real press_rel_error = 0.;
   this->computeContactSolution( nElemsXM, nElemsYM, nElemsZM,
                                 nElemsXS, nElemsYS, nElemsZS,
                                 x_min1, y_min1, z_min1,
                                 x_max1, y_max1, z_max1,
                                 x_min2, y_min2, z_min2,
                                 x_max2, y_max2, z_max2,
                                 thetaM, thetaS,
                                 tribol::SINGLE_MORTAR,
                                 "uniform_mortar_fine",
                                 contact, output, debug, 
                                 visualization, xsol_cntct, press_rel_error );

   this->TearDown();
   this->SetUp();
   contact = false;
   output = false; 
   debug = false;
   visualization = false;

   real press_rel_error_null = 0.;

   this->computeContactSolution( nElemsXM, nElemsYM, nElemsZM, 
                                 nElemsXS, nElemsYS, nElemsZS,
                                 x_min1, y_min1, z_min1,
                                 x_max1, y_max1, z_max1,
                                 x_min2, y_min2, z_min2,
                                 x_max2, y_max2, z_max2,
                                 thetaM, thetaS, 
                                 tribol::SINGLE_MORTAR,
                                 "uniform_mortar_fine_exact",
                                 contact, output, debug, 
                                 visualization, xsol_base, press_rel_error_null );

   mfem::Vector diff(xsol_cntct.Size());
   subtract(xsol_cntct, xsol_base, diff);

   double tol = 5.e-4; // 1/100th of the gap
   int size = this->m_mesh.dim * this->m_mesh.numTotalNodes;
   for (int i=0; i<size; ++i)
   {
      EXPECT_LE( std::abs(diff(i)), tol );
   }

   double press_tol = 1.e-2; // 0.02% error in pressure
   SLIC_INFO( "press_rel_error: " << press_rel_error );
   EXPECT_LE( std::abs(press_rel_error), press_tol );
}

TEST_F( MortarLMPatchTest, single_mortar_nonuniform_mortar_fine_patch )
{

   mfem::Vector xsol_base;
   mfem::Vector xsol_cntct;

   int nMortarElems = 5; // 4; // keep for a nice non-uniform mesh
   int nElemsXM = nMortarElems;
   int nElemsYM = nMortarElems;
   int nElemsZM = 1;

   int nNonmortarElems = 3; // 3; // keep for a nice non-uniform mesh
   int nElemsXS = nNonmortarElems;
   int nElemsYS = nNonmortarElems;
   int nElemsZS = 1;

   real x_min1 = 0.;
   real y_min1 = 0.;
   real z_min1 = 0.; 
   real x_max1 = 1.;
   real y_max1 = 1.;
   real z_max1 = 1.05;

   real x_min2 = 0.;
   real y_min2 = 0.;
   real z_min2 = 0.95;
   real x_max2 = 1.;
   real y_max2 = 1.;
   real z_max2 = 2.;

   bool output = false; // solution vector and rhs
   bool contact = true;
   bool debug = false;  // element matrix writes to .txt files
   bool visualization = false; // output of .vtk with mesh and stress output

   real thetaM = 3.; // 0.; // keep for a nice non-uniform mesh
   real thetaS = -5.; // 3.; // keep for a nice non-uniform mesh

   real press_rel_error = 0.;
   this->computeContactSolution( nElemsXM, nElemsYM, nElemsZM,
                                 nElemsXS, nElemsYS, nElemsZS,
                                 x_min1, y_min1, z_min1,
                                 x_max1, y_max1, z_max1,
                                 x_min2, y_min2, z_min2,
                                 x_max2, y_max2, z_max2,
                                 thetaM, thetaS,
                                 tribol::SINGLE_MORTAR,
                                 "nonuniform_mortar_fine",
                                 contact, output, debug, 
                                 visualization, xsol_cntct, press_rel_error );

   this->TearDown();
   this->SetUp();
   contact = false;
   output = false; 
   debug = false;
   visualization = false;

   real press_rel_error_null = 0.;

   this->computeContactSolution( nElemsXM, nElemsYM, nElemsZM, 
                                 nElemsXS, nElemsYS, nElemsZS,
                                 x_min1, y_min1, z_min1,
                                 x_max1, y_max1, z_max1,
                                 x_min2, y_min2, z_min2,
                                 x_max2, y_max2, z_max2,
                                 thetaM, thetaS, 
                                 tribol::SINGLE_MORTAR,
                                 "nonuniform_mortar_fine_exact",
                                 contact, output, debug, 
                                 visualization, xsol_base, press_rel_error_null );

   mfem::Vector diff(xsol_cntct.Size());
   subtract(xsol_cntct, xsol_base, diff);

   double tol = 5.e-4; // 1/100th of the gap
   int size = this->m_mesh.dim * this->m_mesh.numTotalNodes;
   for (int i=0; i<size; ++i)
   {
      EXPECT_LE( std::abs(diff(i)), tol );
   }

   double press_tol = 5.e-3; // 0.5% error in pressure
   SLIC_INFO( "press_rel_error: " << press_rel_error );
   EXPECT_LE( std::abs(press_rel_error), press_tol );
}

TEST_F( MortarLMPatchTest, single_mortar_nonuniform_nonmortar_fine_patch )
{

   mfem::Vector xsol_base;
   mfem::Vector xsol_cntct;

   int nMortarElems = 3;
   int nElemsXM = nMortarElems;
   int nElemsYM = nMortarElems;
   int nElemsZM = 1;

   int nNonmortarElems = 5;
   int nElemsXS = nNonmortarElems;
   int nElemsYS = 7; //nNonmortarElems;
   int nElemsZS = 1;

   real x_min1 = 0.;
   real y_min1 = 0.;
   real z_min1 = 0.; 
   real x_max1 = 1.;
   real y_max1 = 1.;
   real z_max1 = 1.05;

   real x_min2 = 0.;
   real y_min2 = 0.;
   real z_min2 = 0.95;
   real x_max2 = 1.;
   real y_max2 = 1.;
   real z_max2 = 2.;

   bool output = false; // solution vector and rhs
   bool contact = true;
   bool debug = false;  // element matrix writes to .txt files
   bool visualization = false; // output of .vtk with mesh and stress output

   real press_rel_error = 0.;
   this->computeContactSolution( nElemsXM, nElemsYM, nElemsZM,
                                 nElemsXS, nElemsYS, nElemsZS,
                                 x_min1, y_min1, z_min1,
                                 x_max1, y_max1, z_max1,
                                 x_min2, y_min2, z_min2,
                                 x_max2, y_max2, z_max2,
                                 0., 0.,
                                 tribol::SINGLE_MORTAR,
                                 "uniform_nonmortar_fine",
                                 contact, output, debug, 
                                 visualization, xsol_cntct, press_rel_error );

   this->TearDown();
   this->SetUp();
   contact = false;
   output = false; 
   debug = false;
   visualization = false;

   real press_rel_error_null = 0.;

   this->computeContactSolution( nElemsXM, nElemsYM, nElemsZM, 
                                 nElemsXS, nElemsYS, nElemsZS,
                                 x_min1, y_min1, z_min1,
                                 x_max1, y_max1, z_max1,
                                 x_min2, y_min2, z_min2,
                                 x_max2, y_max2, z_max2,
                                 0., 0.,
                                 tribol::SINGLE_MORTAR,
                                 "uniform_nonmortar_fine_exact",
                                 contact, output, debug, 
                                 visualization, xsol_base, press_rel_error_null );

   mfem::Vector diff(xsol_cntct.Size());
   subtract(xsol_cntct, xsol_base, diff);

   double tol = 1.e-10;
   int size = this->m_mesh.dim * this->m_mesh.numTotalNodes;
   for (int i=0; i<size; ++i)
   {
      EXPECT_LE( std::abs(diff(i)), tol );
   }

   double press_tol = 1.e-7;
   SLIC_INFO( "press_rel_error: " << press_rel_error );
   EXPECT_LE( std::abs(press_rel_error), press_tol );
}

TEST_F( MortarLMPatchTest, aligned_mortar_patch )
{

   mfem::Vector xsol_base;
   mfem::Vector xsol_cntct;

   int nMortarElems = 4; 
   int nElemsXM = nMortarElems;
   int nElemsYM = nMortarElems;
   int nElemsZM = 1;

   int nNonmortarElems = 4;
   int nElemsXS = nNonmortarElems;
   int nElemsYS = nNonmortarElems;
   int nElemsZS = 1;

   real x_min1 = 0.;
   real y_min1 = 0.;
   real z_min1 = 0.; 
   real x_max1 = 1.;
   real y_max1 = 1.;
   real z_max1 = 1.05;

   real x_min2 = 0.;
   real y_min2 = 0.;
   real z_min2 = 0.95;
   real x_max2 = 1.;
   real y_max2 = 1.;
   real z_max2 = 2.;

   bool output = false; // solution output and rhs
   bool contact = true; 
   bool debug = false; // debug matrix writes to .txt files
   bool visualization = false; // visualize mesh and stress

   real press_rel_error = 0.;
   this->computeContactSolution( nElemsXM, nElemsYM, nElemsZM,
                                 nElemsXS, nElemsYS, nElemsZS,
                                 x_min1, y_min1, z_min1,
                                 x_max1, y_max1, z_max1,
                                 x_min2, y_min2, z_min2,
                                 x_max2, y_max2, z_max2,
                                 0., 0., 
                                 tribol::ALIGNED_MORTAR,
                                 "aligned",
                                 contact, output, debug, 
                                 visualization, xsol_cntct, press_rel_error );

   this->TearDown();
   this->SetUp();
   contact = false;
   output = false; 
   debug = false;
   visualization = false;

   real press_rel_error_null = 0.;

   this->computeContactSolution( nElemsXM, nElemsYM, nElemsZM, 
                                 nElemsXS, nElemsYS, nElemsZS,
                                 x_min1, y_min1, z_min1,
                                 x_max1, y_max1, z_max1,
                                 x_min2, y_min2, z_min2,
                                 x_max2, y_max2, z_max2,
                                 0., 0.,
                                 tribol::ALIGNED_MORTAR,
                                 "aligned_exact",
                                 contact, output, debug, 
                                 visualization, xsol_base, press_rel_error_null );

   mfem::Vector diff(xsol_cntct.Size());
   subtract(xsol_cntct, xsol_base, diff);

   double tol = 1.e-10;
   int size = this->m_mesh.dim * this->m_mesh.numTotalNodes;
   for (int i=0; i<size; ++i)
   {
      EXPECT_LE( std::abs(diff(i)), tol );
   }

   double press_tol = 1.e-7;
   SLIC_INFO( "press_rel_error: " << press_rel_error );
   EXPECT_LE( std::abs(press_rel_error), press_tol );
}

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}
