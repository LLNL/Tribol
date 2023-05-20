// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Tribol Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (MIT)

#include "tribol.hpp"

// Tribol includes
#include "tribol/common/logger.hpp"
#include "tribol/common/Parameters.hpp"
#include "tribol/types.hpp"

#include "tribol/mesh/CouplingScheme.hpp"
#include "tribol/mesh/CouplingSchemeManager.hpp"
#include "tribol/mesh/MethodCouplingData.hpp"
#include "tribol/mesh/MeshManager.hpp"
#include "tribol/mesh/MfemData.hpp"
#include "tribol/mesh/InterfacePairs.hpp"

#include "tribol/geom/ContactPlane.hpp"
#include "tribol/geom/ContactPlaneManager.hpp"
#include "tribol/geom/GeomUtilities.hpp"

#include "tribol/physics/Physics.hpp"

#include "tribol/search/InterfacePairFinder.hpp"

#include "tribol/utils/Math.hpp"

// Axom includes
#include "axom/core.hpp"
#include "axom/slic.hpp"

// C/C++ includes
#include <string>
#include <unordered_map>
#include <fstream>

//------------------------------------------------------------------------------
// Interface Implementation
//------------------------------------------------------------------------------
namespace tribol
{

namespace internal
{

/*!
 * \brief Sets various default parameters
 */
void set_defaults()
{
   parameters_t & parameters = parameters_t::getInstance();

   //////////////////////////////
   // visualization parameters //
   //////////////////////////////
   parameters.vis_cycle_incr               = 100;    // default contact output every 100 cycles
   parameters.vis_type                     = VIS_OVERLAPS;
   parameters.output_directory             = "";

   ///////////////////////////////////////
   // computational geometry parameters //
   ///////////////////////////////////////
   parameters.overlap_area_frac            = 1.E-8;  // note: API supports setting this
   parameters.gap_tol_ratio                = 1.e-12; // numerically "zero". Note: API doesn't support setting this
   parameters.gap_separation_ratio         = 0.75;   // applied to the face-radius for geometry filtering. API doesn't support setting this
   parameters.gap_tied_tol                 = 0.1;    // tolerance for how much separation can occur before opposing faces are let go
   parameters.len_collapse_ratio           = 1.E-8;
   parameters.projection_ratio             = 1.E-10;
   parameters.contact_pen_frac             = 3.e-1;  // allows for up to 30% penetration used in timestep vote calculation

}

} /* end namepsace internal */

//------------------------------------------------------------------------------
void initialize( integer dimension, CommType comm )
{
   // sanity checks
   SLIC_ASSERT( (dimension==2) || (dimension==3) );
   SLIC_ASSERT( comm != TRIBOL_COMM_NULL );

   parameters_t & parameters = parameters_t::getInstance();
   parameters.dimension    = dimension;
   parameters.problem_comm = comm;

   internal::set_defaults( );
}

//------------------------------------------------------------------------------
void setPenaltyOptions( int couplingSchemeIndex, PenaltyConstraintType pen_enfrc_option,
                        KinematicPenaltyCalculation kinematic_calc,
                        RatePenaltyCalculation rate_calc )
{
   CouplingSchemeManager& csManager = CouplingSchemeManager::getInstance();

   // check to see if coupling scheme exists
   SLIC_ERROR_IF( !csManager.hasCoupling( couplingSchemeIndex ), 
                  "tribol::setPenaltyOptions(): call tribol::registerCouplingScheme() " <<
                  "prior to calling this routine." );

   // get access to coupling scheme
   CouplingScheme* couplingScheme  = csManager.getCoupling(couplingSchemeIndex);

   // get access to struct on coupling scheme holding penalty options
   EnforcementOptions& enforcement_options = couplingScheme->getEnforcementOptions();
   PenaltyEnforcementOptions& penalty_options = enforcement_options.penalty_options;

   // check that penalty enforcement option is valid
   if ( !in_range(pen_enfrc_option, NUM_PENALTY_OPTIONS) )
   {
      SLIC_WARNING( "tribol::setPenaltyOptions(): penalty enforcement option not available." );
   }
   else
   {
      penalty_options.constraint_type = pen_enfrc_option;
      penalty_options.constraint_type_set = true;
   }

   // check that kinematic penalty calculation is valid
   if ( !in_range(kinematic_calc, NUM_KINEMATIC_PENALTY_CALCULATION) )
   {
      SLIC_WARNING( "tribol::setPenaltyOptions(): kinematic penalty calculation not available." ); 
   }
   else
   {
      penalty_options.kinematic_calculation = kinematic_calc;
      penalty_options.kinematic_calc_set = true;
   }

   // check that the rate penalty calculation is valid
   if ( !in_range(rate_calc, NUM_RATE_PENALTY_CALCULATION) )
   {
      SLIC_WARNING( "tribol::setPenaltyOptions(): rate penalty calculation not available." );
   }
   else
   {
      penalty_options.rate_calculation = rate_calc;
      penalty_options.rate_calc_set = true;
   }

} // end setPenaltyOptions()

//------------------------------------------------------------------------------
void setKinematicConstantPenalty( int meshId, double k )
{
   MeshManager & meshManager = MeshManager::getInstance();

   SLIC_ERROR_IF(!meshManager.hasMesh(meshId), 
                 "tribol::setKinematicConstantPenalty(): " << 
                 "no mesh with id, " << meshId << "exists.");

   registerRealElementField( meshId, KINEMATIC_CONSTANT_STIFFNESS, &k ); 

} // end setKinematicConstantPenalty()

//------------------------------------------------------------------------------
void setKinematicElementPenalty( int meshId, 
                                 double *material_modulus,
                                 double *element_thickness )
{
   MeshManager & meshManager = MeshManager::getInstance();

   SLIC_ERROR_IF(!meshManager.hasMesh(meshId), 
                 "tribol::setKinematicElementPenalty(): " << 
                 "no mesh with id, " << meshId << "exists.");

   SLIC_ERROR_IF(material_modulus == nullptr || element_thickness == nullptr, 
                 "tribol::setKinematicElementPenalty() contains nullptrs for " << 
                 "element_stiffness options.");

   registerRealElementField( meshId, BULK_MODULUS, material_modulus ); 
   registerRealElementField( meshId, ELEMENT_THICKNESS, element_thickness ); 

} // end setKinematicElementPenalty()

//------------------------------------------------------------------------------
void setRateConstantPenalty( int meshId, double r_k )
{
   MeshManager & meshManager = MeshManager::getInstance();

   SLIC_ERROR_IF(!meshManager.hasMesh(meshId), "tribol::setRateConstantPenalty(): " << 
                 "no mesh with id, " << meshId << "exists.");

   registerRealElementField( meshId, RATE_CONSTANT_STIFFNESS, &r_k );

} // end setRateConstantPenalty()

//------------------------------------------------------------------------------
void setRatePercentPenalty( int meshId, double r_p )
{
   MeshManager & meshManager = MeshManager::getInstance();

   SLIC_ERROR_IF(!meshManager.hasMesh(meshId), "tribol::setRatePercentPenalty(): " << 
                 "no mesh with id, " << meshId << "exists.");

   registerRealElementField( meshId, RATE_PERCENT_STIFFNESS, &r_p );

} // end setRatePercentPenalty()

//------------------------------------------------------------------------------
void setContactPenFrac( double frac )
{
   parameters_t & parameters = parameters_t::getInstance();
   if (frac <= 0.)
   {
      // Don't set the contact_pen_frac. This will use default of 30%
      return;
   }

   parameters.contact_pen_frac = frac;

} // end setContactPenFrac()

//------------------------------------------------------------------------------
void setContactAreaFrac( double frac )
{
   parameters_t & parameters = parameters_t::getInstance();
   parameters.overlap_area_frac = frac;
}

//------------------------------------------------------------------------------
void setPenaltyScale( int meshId, double scale )
{
   MeshManager & meshManager = MeshManager::getInstance();

   SLIC_ERROR_IF(!meshManager.hasMesh(meshId), 
                 "tribol::setPenaltyScale(): " << 
                 "no mesh with id, " << meshId << "exists.");

   MeshData & mesh = meshManager.GetMeshInstance( meshId );
   if (scale > 1.e-6)
   {
      mesh.m_elemData.m_penalty_scale = scale;
   }
   else
   {
      // still set small penalty to allow for zeroing out kinematic penalty 
      // enforcement allowing for rate only enforcement
      mesh.m_elemData.m_penalty_scale = scale;
      SLIC_WARNING("tribol::setPenaltyScale(): input scale factor is " << 
                   "close to zero or negative; kinematic contact may " << 
                   "not be properly enforced.");
   }

} // end setPenaltyScale()

//------------------------------------------------------------------------------
void setLagrangeMultiplierOptions( int couplingSchemeIndex, ImplicitEvalMode evalMode, 
                                   SparseMode sparseMode )
{
   // get access to coupling scheme
   CouplingSchemeManager& csManager = CouplingSchemeManager::getInstance();

   SLIC_ERROR_IF( !csManager.hasCoupling( couplingSchemeIndex ), 
                  "tribol::setLagrangeMultiplierOptions(): call tribol::registerCouplingScheme() " <<
                  "prior to calling this routine." );

   CouplingScheme* couplingScheme  = csManager.getCoupling(couplingSchemeIndex);

   // get access to struct on coupling scheme holding penalty options
   EnforcementOptions& enforcement_options = couplingScheme->getEnforcementOptions();
   LagrangeMultiplierImplicitOptions& lm_options = enforcement_options.lm_implicit_options;

   lm_options.eval_mode = evalMode;
   lm_options.sparse_mode = sparseMode;
   if (couplingScheme->getMfemMeshData())
   {
      // MFEM_ELEMENT_DENSE is required to use the MFEM interface
      lm_options.sparse_mode = SparseMode::MFEM_ELEMENT_DENSE;
   }
   lm_options.enforcement_option_set = true;

} // end setLagrangeMultiplierOptions()

//------------------------------------------------------------------------------
void setPlotCycleIncrement( double incr )
{
   parameters_t & parameters = parameters_t::getInstance();
   parameters.vis_cycle_incr = incr;
}

//------------------------------------------------------------------------------
void setPlotOptions( enum VisType v_type )
{
   parameters_t & parameters = parameters_t::getInstance();
   parameters.vis_type = v_type;
}

//------------------------------------------------------------------------------
void setOutputDirectory( const std::string& dir)
{
   // Create path if it doesn't already exist
   if(! axom::utilities::filesystem::pathExists(dir) )
   {
     SLIC_INFO("Creating output path '" << dir << "'");
     axom::utilities::filesystem::makeDirsForPath(dir);
   }

   parameters_t & parameters = parameters_t::getInstance();
   parameters.output_directory = dir;
}

//------------------------------------------------------------------------------
void registerMesh( integer meshId,
                   integer numCells,
                   integer lengthNodalData,
                   const IndexType* connectivity,
                   integer elementType,
                   const real* x,
                   const real* y,
                   const real* z )
{
   MeshManager & meshManager = MeshManager::getInstance();
   MeshData & mesh = meshManager.CreateMesh( meshId );

   // check supported element types
   if (static_cast< InterfaceElementType >(elementType) != EDGE && 
       static_cast< InterfaceElementType >(elementType) != FACE)
   {
      SLIC_WARNING("Mesh topology not supported for mesh id, " << meshId << ".");
      mesh.m_isValid = false;
   }

   const int dim = (z == nullptr) ? 2 : 3;
   if (x == nullptr || y == nullptr)
   {
      SLIC_WARNING("Pointer to x and y-component mesh coordinate arrays are null pointers " <<
                   " for mesh id, " << meshId << ".");
      mesh.m_isValid = false;
   }

   if (dim == 3)
   {
      if (z == nullptr)
      {
         SLIC_WARNING("Pointer to z-component mesh coordinates is null for " << 
                      "mesh id, " << meshId << ".");
         mesh.m_isValid = false;
      }
   }

   // Setup mesh data; input argument pointers are allowed to be null 
   // since Tribol supports null meshes. This is not uncommon in parallel 
   // contact simulations
   mesh.m_meshId = meshId;
   mesh.m_positionX = x;
   mesh.m_positionY = y;
   mesh.m_positionZ = z;
   mesh.m_connectivity = connectivity;
   mesh.m_elementType = static_cast< InterfaceElementType >( elementType );
   mesh.m_lengthNodalData = lengthNodalData;
   mesh.m_numCells = numCells;
  
   // set the number of cells on the mesh element data struct
   mesh.m_elemData.m_numCells = numCells;

   // set the number of nodes on the mesh nodal data struct
   mesh.m_nodalFields.m_numNodes = lengthNodalData;

   // set the number of nodes per cell on the mesh.
   // NOTE: this handles linear segments and bilinear quads
   mesh.m_numCellNodes = (mesh.m_elementType == tribol::EDGE) ? 2 : 4;

   // compute the number of unique surface nodes from the connectivity
   // Note: this routine assigns mesh.m_numSurfaceNodes and allocates
   // space for m_sortedSurfaceNodeIds containing list of unique sorted
   // connectivity node ids in ascending order
   if (mesh.m_numCells > 0)
   {
      mesh.sortSurfaceNodeIds();
   }

   // allocate outward unit face normal arrays and centroid arrays and cell area array
   mesh.m_dim = dim;
   mesh.deallocateArrays();

   if (mesh.m_numCells > 0)
   {
      mesh.allocateArrays(dim);
      initRealArray( mesh.m_nX,   mesh.m_numCells, 0. );
      initRealArray( mesh.m_nY,   mesh.m_numCells, 0. );
      initRealArray( mesh.m_cX,   mesh.m_numCells, 0. );
      initRealArray( mesh.m_cY,   mesh.m_numCells, 0. );
      initRealArray( mesh.m_area, mesh.m_numCells, 0. );
   }

   if (dim == 3 && mesh.m_numCells > 0)
   {
      initRealArray( mesh.m_nZ, mesh.m_numCells, 0. );
      initRealArray( mesh.m_cZ, mesh.m_numCells, 0. );
   }

   // compute the face data
   mesh.computeFaceData( dim );

} // end of registerMesh()

//------------------------------------------------------------------------------
void registerParMesh( integer cs_id,
                      integer mesh_id_1,
                      integer mesh_id_2,
                      mfem::ParMesh& mesh,
                      const mfem::ParGridFunction& current_coords,
                      const std::set<integer>& attributes_1,
                      const std::set<integer>& attributes_2,
                      integer contact_mode,
                      integer contact_case,
                      integer contact_method,
                      integer contact_model,
                      integer enforcement_method,
                      integer binning_method)
{
   std::unique_ptr<mfem::FiniteElementCollection> dual_fec = nullptr;
   integer dual_vdim = 0;
   auto mfem_data = std::make_unique<MfemMeshData>(
      mesh_id_1,
      mesh_id_2,
      mesh,
      current_coords,
      attributes_1,
      attributes_2
   );
   // register empty meshes so the coupling scheme is valid
   registerMesh(
      mesh_id_1, 0, 0, nullptr, mfem_data->GetElemType(), nullptr, nullptr);
   registerMesh(
      mesh_id_2, 0, 0, nullptr, mfem_data->GetElemType(), nullptr, nullptr);
   registerCouplingScheme(
      cs_id,
      mesh_id_1,
      mesh_id_2,
      contact_mode,
      contact_case,
      contact_method,
      contact_model,
      enforcement_method,
      binning_method
   );
   auto coupling_scheme = CouplingSchemeManager::getInstance().getCoupling(cs_id);
   coupling_scheme->setMfemMeshData(std::move(mfem_data));
   if (enforcement_method == LAGRANGE_MULTIPLIER)
   {
      dual_fec = std::make_unique<mfem::H1_FECollection>(
         current_coords.FESpace()->FEColl()->GetOrder(),
         mesh.SpaceDimension()
      );
      if (contact_model == FRICTIONLESS)
      {
         dual_vdim = 1;
      }
      else if (contact_model == TIED || contact_model == COULOMB)
      {
         dual_vdim = mesh.SpaceDimension();
      }
      else
      {
         SLIC_ERROR_ROOT("Unsupported contact model. "
           "Only FRICTIONLESS, TIED, and COULOMB supported.");
      }
      coupling_scheme->setMfemDualData(
         std::make_unique<MfemDualData>(
            mfem_data->GetSubmesh(),
            std::move(dual_fec),
            dual_vdim
         )
      );
   }

} // end of registerParMesh()

//------------------------------------------------------------------------------
void registerNodalDisplacements( integer meshId,
                                 const real* dx,
                                 const real* dy,
                                 const real* dz )
{
   MeshManager & meshManager = MeshManager::getInstance();

   SLIC_ERROR_IF(!meshManager.hasMesh(meshId), "tribol::registerNodalDisplacements(): " << 
                 "no mesh with id, " << meshId << "exists.");

   MeshData & mesh = meshManager.GetMeshInstance( meshId );

   if (dx == nullptr || dy == nullptr)
   {
      mesh.m_isValid = false;
   }

   if (mesh.m_dim == 3)
   {
      if (dz == nullptr)
      {
         mesh.m_isValid = false;
      }
   }

   SLIC_WARNING_IF( !mesh.m_isValid, "tribol::registerNodalDisplacements(): " << 
                    "null pointer to nodal displacement components." );

   mesh.m_dispX = dx;
   mesh.m_dispY = dy;
   mesh.m_dispZ = dz;

} // end registerNodalDisplacements()

//------------------------------------------------------------------------------
void registerNodalVelocities( integer meshId,
                              const real* vx,
                              const real* vy,
                              const real* vz )
{
   MeshManager & meshManager = MeshManager::getInstance();

   SLIC_ERROR_IF(!meshManager.hasMesh(meshId), "tribol::registerNodalVelocities(): " << 
                 "no mesh with id, " << meshId << "exists.");

   MeshData & mesh = meshManager.GetMeshInstance( meshId );
   if (vx == nullptr || vy == nullptr)
   {
      mesh.m_isValid = false;
   }
   
   if (mesh.m_dim == 3)
   {
      if (vz == nullptr)
      {
         mesh.m_isValid = false;
      }
   }   

   SLIC_WARNING_IF( !mesh.m_isValid, "tribol::registerNodalVelocities(): " << 
                    "null pointer to nodal velocity components." );

   mesh.m_velX = vx;
   mesh.m_velY = vy;
   mesh.m_velZ = vz;

   mesh.m_nodalFields.m_is_velocity_set = true;

} // end registerNodalVelocities()

//------------------------------------------------------------------------------
void registerVelocityGridFn( integer cs_id, const mfem::ParGridFunction &v )
{
   CouplingSchemeManager::getInstance().getCoupling(cs_id)->getMfemMeshData()
      ->SetParentVelocity(v);
}

//------------------------------------------------------------------------------
void registerNodalResponse( integer meshId,
                            real* rx,
                            real* ry,
                            real* rz )
{
   MeshManager & meshManager = MeshManager::getInstance();

   SLIC_ERROR_IF(!meshManager.hasMesh(meshId), "tribol::registerNodalResponse(): " << 
                 "no mesh with id, " << meshId << "exists.");

   MeshData & mesh = meshManager.GetMeshInstance( meshId );
   if (rx == nullptr || ry == nullptr)
   {
      mesh.m_isValid = false;
   }

   if (mesh.m_dim == 3)
   {
      if (rz == nullptr)
      {
         mesh.m_isValid = false;
      }
   }   

   SLIC_WARNING_IF(!mesh.m_isValid, "tribol::registerNodalResponse(): " << 
                   "null pointer to nodal response components.");

   mesh.m_forceX = rx;
   mesh.m_forceY = ry;
   mesh.m_forceZ = rz;

} // end registerNodalResponse()

//------------------------------------------------------------------------------
mfem::ParGridFunction getResponseGridFn( integer cs_id )
{
   return CouplingSchemeManager::getInstance().getCoupling(cs_id)
      ->getMfemMeshData()->GetParentResponse();
}

//------------------------------------------------------------------------------
int getMfemSparseMatrix( mfem::SparseMatrix ** sMat, int csId )
{

   if (*sMat != nullptr)
   {
      SLIC_WARNING("tribol::getMfemSparseMatrix(): sparse matrix pointer not null; " <<
                   "nullifying now.");
      *sMat = nullptr;
   }

   CouplingSchemeManager& csManager = CouplingSchemeManager::getInstance();

   SLIC_ERROR_IF(!csManager.hasCoupling(csId), "tribol::getMfemSparseMatrix(): " << 
                 "invalid CouplingScheme id.");

   CouplingScheme* couplingScheme  = csManager.getCoupling( csId );

   switch (couplingScheme->getContactMethod())
   {
      case MORTAR_WEIGHTS:
      case ALIGNED_MORTAR:
      case SINGLE_MORTAR:
      {
         *sMat = static_cast<MortarData*>( couplingScheme->getMethodData() )->getMfemSparseMatrix();
         return 0;
      }
      default:
      {
         SLIC_WARNING("tribol::getMfemSparseMatrix(): interface method does not return matrix data.");
         return 1;
      }
   }
} // end getMfemSparseMatrix()

//------------------------------------------------------------------------------
int getCSRMatrix( int** I, int** J, real** vals, int csId,
                  int* n_offsets, int* n_nonzero )
{
   // check to make sure input pointers are null
   if ( *I != nullptr || *J != nullptr || *vals != nullptr )
   {
      SLIC_WARNING("tribol::getCSRMatrix: input pointers not null, nullifying now.");
      *I = nullptr;
      *J = nullptr;
      *vals = nullptr;
   }

   CouplingSchemeManager& csManager = CouplingSchemeManager::getInstance();

   SLIC_ERROR_IF(!csManager.hasCoupling(csId), "tribol::getCSRMatrix(): invalid " << 
                 "CouplingScheme id.");

   CouplingScheme* couplingScheme  = csManager.getCoupling( csId );

   switch (couplingScheme->getContactMethod())
   {
      case ALIGNED_MORTAR:
      {
         SLIC_WARNING("tribol::getCSRMatrix(): CSR format not currently implemented with " <<
                      "ALIGNED_MORTAR. Use MFEM sparse matrix registration.");
         return 1;
      }
      case MORTAR_WEIGHTS:
      {
         static_cast<MortarData*>( couplingScheme->getMethodData() )->getCSRArrays( I, J, vals, n_offsets, n_nonzero );
         return 0;
      }
      case SINGLE_MORTAR:
      {
         SLIC_WARNING("tribol::getCSRMatrix(): CSR format not currently implemented with "
                      "SINGLE_MORTAR. Use MFEM sparse matrix registration.");
         return 1;
      }
      default:
      {
         SLIC_WARNING("tribol::registerCSRMatrix(): method does not return matrix data; " <<
                       "invalid call.");
         return 1;
      }
   }
} // end getCSRMatrix()

//------------------------------------------------------------------------------
int getElementBlockJacobians( integer csId, 
                              BlockSpace row_block, 
                              BlockSpace col_block,
                              const axom::Array<integer>* row_elem_idx,
                              const axom::Array<integer>* col_elem_idx,
                              const axom::Array<mfem::DenseMatrix>* jacobians )
{
   SparseMode sparse_mode = CouplingSchemeManager::getInstance().
      getCoupling(csId)->getEnforcementOptions().lm_implicit_options.sparse_mode;
   if (sparse_mode != SparseMode::MFEM_ELEMENT_DENSE)
   {
      SLIC_WARNING("Jacobian is assembled and can be accessed by " 
         "getMfemSparseMatrix() or getCSRMatrix(). For (unassembled) element "
         "Jacobian contributions, call setLagrangeMultiplierOptions() with "
         "SparseMode::MFEM_ELEMENT_DENSE before calling update().");
      return 1;
   }
   MethodData* method_data = 
      CouplingSchemeManager::getInstance().getCoupling( csId )->getMethodData();
   row_elem_idx = &method_data->getBlockJElementIds()[static_cast<int>(row_block)];
   col_elem_idx = &method_data->getBlockJElementIds()[static_cast<int>(col_block)];
   jacobians = &method_data->getBlockJ()(
      static_cast<int>(row_block),
      static_cast<int>(col_block)
   );
   return 0;
}

//------------------------------------------------------------------------------
std::unique_ptr<mfem::BlockOperator> getMfemBlockJacobian( integer csId )
{
   CouplingScheme* coupling_scheme = CouplingSchemeManager::getInstance().
      getCoupling(csId);
   SparseMode sparse_mode = coupling_scheme
      ->getEnforcementOptions().lm_implicit_options.sparse_mode;
   if (sparse_mode != SparseMode::MFEM_ELEMENT_DENSE)
   {
      SLIC_ERROR_ROOT("Jacobian is assembled and can be accessed by " 
         "getMfemSparseMatrix() or getCSRMatrix(). For (unassembled) element "
         "Jacobian contributions, call setLagrangeMultiplierOptions() with "
         "SparseMode::MFEM_ELEMENT_DENSE before calling update().");
   }
   auto mfem_data = coupling_scheme->getMfemMeshData();
   SLIC_ERROR_ROOT_IF(!mfem_data, "No MFEM data exists. The coupling scheme "
     "must be registered using registerParMesh() to use this method.");
   MethodData* method_data = coupling_scheme->getMethodData();
   // 0 = displacement DOFs, 1 = lagrange multiplier DOFs
   // (0,0) block is empty (for now)
   // (1,1) block is empty
   // (0,1) and (1,0) are symmetric (for now)
   const auto& elem_map_1 = mfem_data->GetElemMap1();
   const auto& elem_map_2 = mfem_data->GetElemMap2();
   auto mortar_elems = method_data
      ->getBlockJElementIds()[static_cast<int>(BlockSpace::MORTAR)];
   for (auto& mortar_elem : mortar_elems)
   {
      mortar_elem = elem_map_1[static_cast<size_t>(mortar_elem)];
   }
   auto nonmortar_elems = method_data
      ->getBlockJElementIds()[static_cast<int>(BlockSpace::NONMORTAR)];
   for (auto& nonmortar_elem : nonmortar_elems)
   {
      nonmortar_elem = elem_map_2[static_cast<size_t>(nonmortar_elem)];
   }
   auto lm_elems = method_data
      ->getBlockJElementIds()[static_cast<int>(BlockSpace::LAGRANGE_MULTIPLIER)];
   for (auto& lm_elem : lm_elems)
   {
      lm_elem = elem_map_2[static_cast<size_t>(lm_elem)];
   }
   // get (1,0) block
   const auto& elem_J_1 = method_data->getBlockJ()(
      static_cast<int>(BlockSpace::LAGRANGE_MULTIPLIER),
      static_cast<int>(BlockSpace::MORTAR)
   );
   const auto& elem_J_2 = method_data->getBlockJ()(
      static_cast<int>(BlockSpace::LAGRANGE_MULTIPLIER),
      static_cast<int>(BlockSpace::NONMORTAR)
   );
   // move to submesh level
   auto matrix_xfer = coupling_scheme->getMatrixXfer();
   auto submesh_J = matrix_xfer->TransferToParallelSparse(
      lm_elems, 
      mortar_elems, 
      elem_J_1
   );
   submesh_J += matrix_xfer->TransferToParallelSparse(
      lm_elems, 
      nonmortar_elems, 
      elem_J_2
   );
   submesh_J.Finalize();

   // transform J values from submesh to parent
   auto J = submesh_J.GetJ();
   for (int j{}; j < submesh_J.NumNonZeroElems(); ++j)
   {
      J[j] = (*coupling_scheme->getSubmeshToParentVdofList())[J[j]];
   }

   // create block operator
   auto block_J = std::make_unique<mfem::BlockOperator>(
      *coupling_scheme->getBlockOffsets()
   );
   block_J->owns_blocks = 1;

   // fill block operator
   auto hypre_J = matrix_xfer->ConvertToHypreParMatrix(submesh_J);
   block_J->SetBlock(0, 1, new mfem::TransposeOperator(hypre_J.get()));
   block_J->SetBlock(1, 0, hypre_J.release());

   return block_J;
}

//------------------------------------------------------------------------------
void registerMortarGaps( integer meshId,
                         real * nodal_gaps )
{
   MeshManager & meshManager = MeshManager::getInstance();

   SLIC_ERROR_IF(!meshManager.hasMesh(meshId), "tribol::registerMortarGaps(): " << 
                 "no mesh with id " << meshId << " exists.");

   MeshData & mesh = meshManager.GetMeshInstance( meshId );

   if (nodal_gaps == nullptr)
   {
      SLIC_WARNING( "tribol::registerMortarGaps(): null pointer to data " << 
                    "on mesh " << meshId << ".");
      mesh.m_isValid = false;
   }
   else
   {
      mesh.m_nodalFields.m_node_gap = nodal_gaps;
      mesh.m_nodalFields.m_is_node_gap_set = true;
   }
   
}

//------------------------------------------------------------------------------
mfem::ParGridFunction getGapGridFn( integer cs_id )
{
   return CouplingSchemeManager::getInstance().getCoupling(cs_id)
      ->getMfemDualData()->GetSubmeshGap();
}

//------------------------------------------------------------------------------
void registerMortarPressures( integer meshId,
                              const real * nodal_pressures )
{
   MeshManager & meshManager = MeshManager::getInstance();

   SLIC_ERROR_IF(!meshManager.hasMesh(meshId), "tribol::registerMortarPressures(): " << 
                 "no mesh with id " << meshId << " exists.");

   MeshData & mesh = meshManager.GetMeshInstance( meshId );

   if (nodal_pressures == nullptr)
   {
      SLIC_WARNING( "tribol::registerMortarPressures(): null pointer to data " << 
                    "on mesh " << meshId << ".");
      mesh.m_isValid = false;
   }
   else
   {
      mesh.m_nodalFields.m_node_pressure = nodal_pressures;
      mesh.m_nodalFields.m_is_node_pressure_set = true;
   }
   
}

//------------------------------------------------------------------------------
mfem::ParGridFunction& getPressureGridFn( integer cs_id )
{
   return CouplingSchemeManager::getInstance().getCoupling(cs_id)
      ->getMfemDualData()->GetSubmeshPressure();
}

//------------------------------------------------------------------------------
void registerIntNodalField( integer meshId,
                            const IntNodalFields field,
                            integer * TRIBOL_UNUSED_PARAM(fieldVariable) )
{
   MeshManager & meshManager = MeshManager::getInstance();

   SLIC_ERROR_IF(!meshManager.hasMesh(meshId), "tribol::registerIntNodalField(): " << 
                 "no mesh with id, " << meshId << "exists.");

   switch (field)
   {
      case UNDEFINED_INT_NODAL_FIELD:
      default:
         SLIC_WARNING("tribol::registerIntNodalField() not yet implemented.");
   } // end switch over field

} // end registerIntNodalField()

//------------------------------------------------------------------------------
void registerRealElementField( integer meshId,
                               const RealElementFields field,
                               real * fieldVariable )
{
   MeshManager & meshManager = MeshManager::getInstance();

   SLIC_ERROR_IF(!meshManager.hasMesh(meshId), "tribol::registerRealElementField(): " << 
               "no mesh with id, " << meshId << "exists.");

   MeshData & mesh = meshManager.GetMeshInstance( meshId );

   switch (field)
   {
      case KINEMATIC_CONSTANT_STIFFNESS:
      {
         if (fieldVariable==nullptr)
         {
            SLIC_WARNING( "tribol::registerRealElementField(): null pointer to data for " << 
                          "'KINEMATIC_CONSTANT_STIFFNESS' on mesh, " << meshId << ".");
            mesh.m_isValid = false;
         }
         else
         {
            mesh.m_elemData.m_penalty_stiffness = *fieldVariable;
            mesh.m_elemData.m_is_kinematic_constant_penalty_set = true;
         }
         break;
      }
      case RATE_CONSTANT_STIFFNESS:
      {
         if (fieldVariable==nullptr)
         {
            SLIC_WARNING( "tribol::registerRealElementField(): null pointer to data for " << 
                          "'RATE_CONSTANT_STIFFNESS' on mesh, " << meshId << ".");
            mesh.m_isValid = false;
         }
         else
         {
            mesh.m_elemData.m_rate_penalty_stiffness = *fieldVariable;
            mesh.m_elemData.m_is_rate_constant_penalty_set = true;
         }
         break;
      }
      case RATE_PERCENT_STIFFNESS:
      {
         if (fieldVariable==nullptr)
         {
            SLIC_WARNING( "tribol::registerRealElementField(): null pointer to data for " << 
                          "'RATE_PERCENT_STIFFNESS' on mesh, " << meshId << ".");
            mesh.m_isValid = false;
         }
         else
         {
            mesh.m_elemData.m_rate_percent_stiffness = *fieldVariable;
            mesh.m_elemData.m_is_rate_percent_penalty_set = true;
         }
         break;
      }
      case BULK_MODULUS:
      {
         if (fieldVariable==nullptr)
         {
            SLIC_WARNING( "tribol::registerRealElementField(): null pointer to data for " << 
                          "'BULK_MODULUS' on mesh, " << meshId << ".");
            mesh.m_isValid = false;
         }
         else
         {
            mesh.m_elemData.m_mat_mod = fieldVariable;
         }

         if (mesh.m_elemData.m_thickness != nullptr)
         {
            mesh.m_elemData.m_is_kinematic_element_penalty_set = true;
         }
         break;
      }
      case YOUNGS_MODULUS:
      {
         if (fieldVariable==nullptr)
         {
            SLIC_WARNING( "tribol::registerRealElementField(): null pointer to data for " << 
                          "'YOUNGS_MODULUS' on mesh, " << meshId << ".");
            mesh.m_isValid = false;
         }
         else
         {
            mesh.m_elemData.m_mat_mod = fieldVariable;
         }

         if (mesh.m_elemData.m_thickness != nullptr)
         {
            mesh.m_elemData.m_is_kinematic_element_penalty_set = true;
         }
         break;
      }
      case ELEMENT_THICKNESS:
      {
         if (fieldVariable==nullptr)
         {
            SLIC_WARNING( "tribol::registerRealElementField(): null pointer to data for " << 
                          "'ELEMENT_THICKNESS' on mesh, " << meshId << ".");
            mesh.m_isValid = false;
         }
         else
         {
            mesh.m_elemData.m_thickness = fieldVariable;
         }

         if (mesh.m_elemData.m_mat_mod != nullptr)
         {
            mesh.m_elemData.m_is_kinematic_element_penalty_set = true;
         }
         break;
      }
      default:
      {
         SLIC_WARNING( "tribol::registerRealElementField(): the field argument " << 
                       "on mesh," << meshId << ", is not an accepted tribol real element field." );
      }
   } // end switch over field

} // end registerRealElementField()

//------------------------------------------------------------------------------
void registerIntElementField( integer meshId,
                              const IntElementFields field,
                              integer * TRIBOL_UNUSED_PARAM(fieldVariable) )
{
   MeshManager & meshManager = MeshManager::getInstance();

   SLIC_ERROR_IF(!meshManager.hasMesh(meshId), "tribol::registerIntElementField(): " << 
                 "no mesh with id, " << meshId << "exists.");

   switch (field)
   {
      case UNDEFINED_INT_ELEMENT_FIELD:
      default:
         SLIC_WARNING("tribol::registerIntElementField() not yet implemented.");
   } // end switch over field

} // end registerIntElementField()

//------------------------------------------------------------------------------
void registerCouplingScheme( integer couplingSchemeIndex,
                             integer meshId1,
                             integer meshId2,
                             integer contact_mode,
                             integer contact_case,
                             integer contact_method,
                             integer contact_model,
                             integer enforcement_method,
                             integer binning_method )
{
   // add coupling scheme. Checks for valid schemes will be performed later
   CouplingSchemeManager& couplingSchemeManager =
         CouplingSchemeManager::getInstance();

   CouplingScheme* scheme =
         new CouplingScheme(couplingSchemeIndex,
                            meshId1,
                            meshId2,
                            contact_mode,
                            contact_case,
                            contact_method,
                            contact_model,
                            enforcement_method,
                            binning_method);

   // add coupling scheme to manager. Validity checks are performed in 
   // tribol::update() when each coupling scheme is initialized.
   couplingSchemeManager.addCoupling(couplingSchemeIndex, scheme);

} // end registerCouplingScheme()

//------------------------------------------------------------------------------
void setInterfacePairs( integer couplingSchemeIndex,
                        IndexType numPairs,
                        IndexType const * const meshId1,
                        IndexType const * const pairType1,
                        IndexType const * const pairIndex1,
                        IndexType const * const meshId2,
                        IndexType const * const pairType2,
                        IndexType const * const pairIndex2 )
{
   CouplingSchemeManager& csManager = CouplingSchemeManager::getInstance();

   SLIC_ERROR_IF(!csManager.hasCoupling(couplingSchemeIndex), 
                 "tribol::setInterfacePairs(): invalid coupling scheme index.");

   auto* couplingScheme = csManager.getCoupling(couplingSchemeIndex);
   auto* pairs = couplingScheme->getInterfacePairs();

   pairs->clear();

   // copy the interaction pairs
   for(int i=0; i< numPairs; ++i)
   {
      // Perform basic geometry proximity filter prior to adding interface pair
      // since a user may specify interface pairs that aren't, nor should be in contact.
      InterfacePair pair { meshId1[i], pairType1[i], pairIndex1[i],
                           meshId2[i], pairType2[i], pairIndex2[i],
                           false, i };
      ContactMode mode = couplingScheme->getContactMode();
      bool check = geomFilter( pair, mode );

      if (check)
      {
         pairs->addInterfacePair( pair );
      }
   }

   // Disable per-cycle rebinning
   couplingScheme->setFixedBinning(true);
}

//------------------------------------------------------------------------------
integer update( integer cycle, real t, real &dt )
{
   CouplingSchemeManager& csManager = CouplingSchemeManager::getInstance();
   int numCouplings = csManager.getNumberOfCouplings();
   bool err_cs = false;

   /////////////////////////////////////////////////////////////////////////
   //                                                                     //
   // Loop over coupling schemes.                                         //
   //                                                                     //
   // Note, numCouplings is always 1 greater than the highest registered  //
   // coupling index. This allows for non-contiguous coupling scheme ids, //
   // which may arise from host-code registration or from skipped schemes // 
   //                                                                     //
   /////////////////////////////////////////////////////////////////////////
   for(int csIndex =0; csIndex < numCouplings; ++csIndex)
   {
      if(!csManager.hasCoupling(csIndex))
      {
         continue;
      }

      CouplingScheme* couplingScheme  = csManager.getCoupling(csIndex);

      // update redecomp meshes if supplied mfem data
      if (couplingScheme->hasMfemData())
      {
         auto mfem_data = couplingScheme->getMfemMeshData();
         auto mesh_id_1 = mfem_data->GetMesh1ID();
         auto mesh_id_2 = mfem_data->GetMesh2ID();
         mfem_data->UpdateMeshData();
         auto coord_ptrs = mfem_data->GetRedecompCoordsPtrs();
         registerMesh(
            mesh_id_1,
            mfem_data->GetMesh1NE(),
            mfem_data->GetNV(),
            mfem_data->GetMesh1Conn(),
            mfem_data->GetElemType(),
            coord_ptrs[0],
            coord_ptrs[1],
            coord_ptrs[2]
         );
         registerMesh(
            mesh_id_2,
            mfem_data->GetMesh2NE(),
            mfem_data->GetNV(),
            mfem_data->GetMesh2Conn(),
            mfem_data->GetElemType(),
            coord_ptrs[0],
            coord_ptrs[1],
            coord_ptrs[2]
         );
         auto f_ptrs = mfem_data->GetRedecompResponsePtrs();
         registerNodalResponse(
            mesh_id_1, f_ptrs[0], f_ptrs[1], f_ptrs[2]);
         registerNodalResponse(
            mesh_id_2, f_ptrs[0], f_ptrs[1], f_ptrs[2]);
         if (mfem_data->HasVelocity())
         {
            auto v_ptrs = mfem_data->GetRedecompVelocityPtrs();
            registerNodalVelocities(
               mesh_id_1, v_ptrs[0], v_ptrs[1], v_ptrs[2]);
            registerNodalVelocities(
               mesh_id_2, v_ptrs[0], v_ptrs[1], v_ptrs[2]);
         }
         if (couplingScheme->getEnforcementMethod() == LAGRANGE_MULTIPLIER)
         {
            SLIC_ERROR_ROOT_IF(couplingScheme->getContactModel() != FRICTIONLESS,
              "Only frictionless contact is supported at this time.");
            auto dual_data = couplingScheme->getMfemDualData();
            dual_data->UpdateDualData(mfem_data->GetRedecompMesh());
            auto g_ptrs = dual_data->GetRedecompGapPtrs();
            registerMortarGaps(mesh_id_2, g_ptrs[0]);
            auto p_ptrs = dual_data->GetRedecompPressurePtrs();
            registerMortarPressures(mesh_id_2, p_ptrs[0]);
            auto lm_options = couplingScheme->getEnforcementOptions().lm_implicit_options;
            if (
               lm_options.enforcement_option_set && 
               (
                  lm_options.eval_mode == ImplicitEvalMode::MORTAR_JACOBIAN ||
                  lm_options.eval_mode == ImplicitEvalMode::MORTAR_RESIDUAL_JACOBIAN
               )
            )
            {
               couplingScheme->setMatrixXfer();
            }
         }
      }
      
      // initialize and check for valid coupling scheme
      if (!couplingScheme->init())
      {
         // should we error out entirely for one invalid coupling scheme? SRW
         SLIC_WARNING("Invalid coupling scheme; please see warnings.");
         continue;
      }

      // perform binning between meshes on the coupling scheme
      couplingScheme->performBinning();

      // check to make sure binning returns non-null pairs
      auto* pairs = couplingScheme->getInterfacePairs();
      if(pairs == nullptr)
      {
         continue;
      }

      err_cs = couplingScheme->apply( cycle, t, dt );

      if ( err_cs != 0 )
      {
         SLIC_WARNING("tribol::update(): coupling scheme, " << csIndex <<
                      ", returned with an error.");
      }

   } // end of coupling scheme loop

   return err_cs;

} // end update()

//------------------------------------------------------------------------------
void finalize( )
{
   CouplingSchemeManager& csManager = CouplingSchemeManager::getInstance();
   int numCouplings = csManager.getNumberOfCouplings();
   if (numCouplings > 0)
   {
      csManager.clearAllCouplings();
   }
}

//------------------------------------------------------------------------------
} // end tribol namespace

