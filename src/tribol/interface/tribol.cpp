// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Tribol Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (MIT)

#include "tribol.hpp"

// Tribol includes
#include "tribol/common/Parameters.hpp"
#include "tribol/types.hpp"

#include "tribol/mesh/CouplingScheme.hpp"
#include "tribol/mesh/MethodCouplingData.hpp"
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
   parameters.auto_contact_pen_frac        = 0.95;   // max allowable interpenetration as percent of element thickness for contact candidacy 
   parameters.timestep_pen_frac            = 3.e-1;  // max allowable interpenetration as percent of element thickness prior to triggering timestep vote (not exposed to API) 
   parameters.enable_timestep_vote         = false;  // true if host-code wants to receive tribol timestep vote
   
   // Interpenetration check for auto-contact. If true, this will check a full-overlap 
   // face-pair configuration in the computational geoemtry routines to preclude 
   // auto-contact of opposite sides of thin structures/plates. If the full-overlap 
   // interpenetration kinematic gap is more than the smallest thickness of the 
   // constituent face elements, then we don't consider the face-pair a contact candidate.
   // Note, auto-contact will require registration of element thicknesses.
   parameters.auto_interpen_check           = false; // true if the auto-contact interpenetration check is used for interpenetrating face-pairs.

}

} /* end namepsace internal */

//------------------------------------------------------------------------------
void initialize( int dimension, CommT comm )
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
void setPenaltyOptions( IndexT cs_id, PenaltyConstraintType pen_enfrc_option,
                        KinematicPenaltyCalculation kinematic_calc,
                        RatePenaltyCalculation rate_calc )
{
   auto couplingScheme = CouplingSchemeManager::getInstance().findData(cs_id);

   // check to see if coupling scheme exists
   SLIC_ERROR_ROOT_IF( !couplingScheme, 
                       "tribol::setPenaltyOptions(): call tribol::registerCouplingScheme() " <<
                       "prior to calling this routine." );

   // get access to struct on coupling scheme holding penalty options
   EnforcementOptions& enforcement_options = couplingScheme->getEnforcementOptions();
   PenaltyEnforcementOptions& penalty_options = enforcement_options.penalty_options;

   // check that penalty enforcement option is valid
   if ( !in_range(pen_enfrc_option, NUM_PENALTY_OPTIONS) )
   {
      SLIC_WARNING_ROOT( "tribol::setPenaltyOptions(): penalty enforcement option not available." );
   }
   else
   {
      penalty_options.constraint_type = pen_enfrc_option;
      penalty_options.constraint_type_set = true;
   }

   // check that kinematic penalty calculation is valid
   if ( !in_range(kinematic_calc, NUM_KINEMATIC_PENALTY_CALCULATION) )
   {
      SLIC_WARNING_ROOT( "tribol::setPenaltyOptions(): kinematic penalty calculation not available." ); 
   }
   else
   {
      penalty_options.kinematic_calculation = kinematic_calc;
      penalty_options.kinematic_calc_set = true;
   }

   // check that the rate penalty calculation is valid
   if ( !in_range(rate_calc, NUM_RATE_PENALTY_CALCULATION) )
   {
      SLIC_WARNING_ROOT( "tribol::setPenaltyOptions(): rate penalty calculation not available." );
   }
   else
   {
      penalty_options.rate_calculation = rate_calc;
      penalty_options.rate_calc_set = true;
   }

} // end setPenaltyOptions()

//------------------------------------------------------------------------------
void setKinematicConstantPenalty( IndexT mesh_id, RealT k )
{
   // note, error checking done in the following registration routine
   registerRealElementField( mesh_id, KINEMATIC_CONSTANT_STIFFNESS, &k ); 

} // end setKinematicConstantPenalty()

//------------------------------------------------------------------------------
void setKinematicElementPenalty( IndexT mesh_id, 
                                 const RealT *material_modulus,
                                 const RealT *element_thickness )
{
   // note, error checking done in the following registration routine
   registerRealElementField( mesh_id, BULK_MODULUS, material_modulus ); 
   registerRealElementField( mesh_id, ELEMENT_THICKNESS, element_thickness ); 

} // end setKinematicElementPenalty()

//------------------------------------------------------------------------------
void setRateConstantPenalty( IndexT mesh_id, RealT r_k )
{
   // note, error checking done in the following registration routine
   registerRealElementField( mesh_id, RATE_CONSTANT_STIFFNESS, &r_k );

} // end setRateConstantPenalty()

//------------------------------------------------------------------------------
void setRatePercentPenalty( IndexT mesh_id, RealT r_p )
{
   // note, error checking done in the following registration routine
   registerRealElementField( mesh_id, RATE_PERCENT_STIFFNESS, &r_p );

} // end setRatePercentPenalty()

//------------------------------------------------------------------------------
void setAutoContactPenScale( RealT scale )
{
   parameters_t & parameters = parameters_t::getInstance();

   // check for strict positivity of the input parameter
   SLIC_WARNING_ROOT_IF(scale<0., "tribol::setAutoContactPenScale(): " << 
                        "input for the auto-contact length scale factor must be positive.");

   parameters.auto_contact_pen_frac = scale;

} // end setAutoContactPenScale()

//------------------------------------------------------------------------------
void setTimestepPenFrac( RealT frac )
{
   parameters_t & parameters = parameters_t::getInstance();
   if (frac <= 0.)
   {
      // Don't set the timestep_pen_frac. This will use default
      return;
   }

   parameters.timestep_pen_frac = frac;

} // end setTimestepPenFrac()

//------------------------------------------------------------------------------
void setContactAreaFrac( RealT frac )
{
   parameters_t & parameters = parameters_t::getInstance();
   if (frac < 1.e-12)
   {
      SLIC_DEBUG_ROOT("tribol::setContactAreaFrac(): area fraction too small or negative; " << 
                      "setting to default 1.e-8.");
      frac = 1.e-8;
   }
   parameters.overlap_area_frac = frac;

} // end setPenaltyScale()

//------------------------------------------------------------------------------
void setPenaltyScale( IndexT mesh_id, RealT scale )
{
   auto mesh = MeshManager::getInstance().findData(mesh_id);

   SLIC_ERROR_ROOT_IF(!mesh, 
                      "tribol::setPenaltyScale(): " << 
                      "no mesh with id, " << mesh_id << "exists.");

   if (scale > 1.e-6)
   {
      mesh->getElementData().m_penalty_scale = scale;
   }
   else
   {
      // still set small penalty to allow for zeroing out kinematic penalty 
      // enforcement allowing for rate only enforcement
      mesh->getElementData().m_penalty_scale = scale;
      SLIC_WARNING_ROOT("tribol::setPenaltyScale(): input scale factor is " << 
                        "close to zero or negative; kinematic contact may " << 
                        "not be properly enforced.");
   }

} // end setPenaltyScale()

//------------------------------------------------------------------------------
void setLagrangeMultiplierOptions( IndexT cs_id, ImplicitEvalMode evalMode, 
                                   SparseMode sparseMode )
{
   // get access to coupling scheme
   auto couplingScheme = CouplingSchemeManager::getInstance().findData(cs_id);

   SLIC_ERROR_ROOT_IF( !couplingScheme, 
                       "tribol::setLagrangeMultiplierOptions(): call tribol::registerCouplingScheme() " <<
                       "prior to calling this routine." );

   // get access to struct on coupling scheme holding penalty options
   EnforcementOptions& enforcement_options = couplingScheme->getEnforcementOptions();
   LagrangeMultiplierImplicitOptions& lm_options = enforcement_options.lm_implicit_options;

   lm_options.eval_mode = evalMode;
   lm_options.sparse_mode = sparseMode;

#ifdef BUILD_REDECOMP

   if (couplingScheme->hasMfemData())
   {
      // MFEM_ELEMENT_DENSE is required to use the MFEM interface
      lm_options.sparse_mode = SparseMode::MFEM_ELEMENT_DENSE;
      if (
         !couplingScheme->hasMfemJacobianData() && (
            lm_options.eval_mode == ImplicitEvalMode::MORTAR_JACOBIAN ||
            lm_options.eval_mode == ImplicitEvalMode::MORTAR_RESIDUAL_JACOBIAN
         )
      )
      {
         couplingScheme->setMfemJacobianData(std::make_unique<MfemJacobianData>(
            *couplingScheme->getMfemMeshData(),
            *couplingScheme->getMfemSubmeshData(),
            couplingScheme->getContactMethod()
         ));
      }
   }

#endif /* BUILD_REDECOMP */

   lm_options.enforcement_option_set = true;

} // end setLagrangeMultiplierOptions()

//------------------------------------------------------------------------------
void setPlotCycleIncrement( int incr )
{
   parameters_t & parameters = parameters_t::getInstance();
   parameters.vis_cycle_incr = incr;

} // end setPlotCycleIncrement()

//------------------------------------------------------------------------------
void setPlotOptions( enum VisType v_type )
{
   parameters_t & parameters = parameters_t::getInstance();
   parameters.vis_type = v_type;

} // end setPlotOptions()

//------------------------------------------------------------------------------
void setOutputDirectory( const std::string& dir)
{
   // Create path if it doesn't already exist
   if(! axom::utilities::filesystem::pathExists(dir) )
   {
     SLIC_INFO_ROOT("Creating output path '" << dir << "'");
     axom::utilities::filesystem::makeDirsForPath(dir);
   }

   parameters_t & parameters = parameters_t::getInstance();
   parameters.output_directory = dir;

} // end setOutputDirectory()

//------------------------------------------------------------------------------
void setLoggingLevel( IndexT cs_id, LoggingLevel log_level )
{
   // get access to coupling scheme
   auto couplingScheme = CouplingSchemeManager::getInstance().findData(cs_id);

   SLIC_ERROR_IF(!couplingScheme, "tribol::setLoggingLevel(): " << 
                 "invalid CouplingScheme id.");

   if ( !in_range(static_cast<int>(log_level), 
                  static_cast<int>(tribol::NUM_LOGGING_LEVELS)) )
   {
      SLIC_INFO_ROOT("tribol::setLoggingLevel(): Logging level not an option; " << 
                     "using 'warning' level.");
      couplingScheme->setLoggingLevel( tribol::TRIBOL_WARNING );
   }
   else
   {
      couplingScheme->setLoggingLevel( log_level );
   }

} // end setLoggingLevel()

//------------------------------------------------------------------------------
void enableTimestepVote( const bool enable )
{
   parameters_t & parameters = parameters_t::getInstance();
   parameters.enable_timestep_vote= enable;

} // end enableTimestepVote()

//------------------------------------------------------------------------------
void registerMesh( IndexT mesh_id,
                   IndexT num_cells,
                   IndexT num_nodes,
                   const IndexT* connectivity,
                   int element_type,
                   const RealT* x,
                   const RealT* y,
                   const RealT* z )
{
   MeshManager::getInstance().addData(mesh_id, MeshData(
      mesh_id, num_cells, num_nodes, connectivity, 
      static_cast<InterfaceElementType>(element_type), x, y, z));
} // end registerMesh()

//------------------------------------------------------------------------------
void registerNodalDisplacements( IndexT mesh_id,
                                 const RealT* dx,
                                 const RealT* dy,
                                 const RealT* dz )
{
   auto mesh = MeshManager::getInstance().findData(mesh_id);

   SLIC_ERROR_ROOT_IF(!mesh, "tribol::registerNodalDisplacements(): " << 
                      "no mesh with id, " << mesh_id << "exists.");

   mesh->getNodalFields().m_is_nodal_displacement_set = true;

   if (dx == nullptr || dy == nullptr)
   {
      mesh->getNodalFields().m_is_nodal_displacement_set = false;
   }

   if (mesh->dimension() == 3)
   {
      if (dz == nullptr)
      {
         mesh->getNodalFields().m_is_nodal_displacement_set = false;
      }
   }

   mesh->setDisplacement(mesh->numberOfNodes(), dx, dy, dz);

} // end registerNodalDisplacements()

//------------------------------------------------------------------------------
void registerNodalVelocities( IndexT mesh_id,
                              const RealT* vx,
                              const RealT* vy,
                              const RealT* vz )
{
   auto mesh = MeshManager::getInstance().findData(mesh_id);

   SLIC_ERROR_ROOT_IF(!mesh, "tribol::registerNodalVelocities(): " << 
                      "no mesh with id, " << mesh_id << "exists.");

   mesh->getNodalFields().m_is_velocity_set = true;

   if (vx == nullptr || vy == nullptr)
   {
      mesh->getNodalFields().m_is_velocity_set = false;
   }
   
   if (mesh->dimension() == 3)
   {
      if (vz == nullptr)
      {
         mesh->getNodalFields().m_is_velocity_set = false;
      }
   }

   mesh->setVelocity(mesh->numberOfNodes(), vx, vy, vz);

} // end registerNodalVelocities()

//------------------------------------------------------------------------------
void registerNodalResponse( IndexT mesh_id,
                            RealT* rx,
                            RealT* ry,
                            RealT* rz )
{
   auto mesh = MeshManager::getInstance().findData(mesh_id);

   SLIC_ERROR_ROOT_IF(!mesh, "tribol::registerNodalResponse(): " << 
                      "no mesh with id, " << mesh_id << "exists.");

   mesh->getNodalFields().m_is_nodal_response_set = true;
   if (rx == nullptr || ry == nullptr)
   {
      mesh->getNodalFields().m_is_nodal_response_set = false;
   }

   if (mesh->dimension() == 3)
   {
      if (rz == nullptr)
      {
         mesh->getNodalFields().m_is_nodal_response_set = false;
      }
   }   

   mesh->setResponse(mesh->numberOfNodes(), rx, ry, rz);

} // end registerNodalResponse()

//------------------------------------------------------------------------------
int getJacobianSparseMatrix( mfem::SparseMatrix ** sMat, IndexT cs_id )
{
   // note, SLIC_ERROR_ROOT_IF is not used here because it's possible not all ranks 
   // will have method (i.e. mortar) data.
   SLIC_ERROR_IF(*sMat!=nullptr, "tribol::getJacobianSparseMatrix(): " << 
                 "sparse matrix pointer not null.");

   // get access to coupling scheme
   auto couplingScheme = CouplingSchemeManager::getInstance().findData(cs_id);

   SLIC_ERROR_IF(!couplingScheme, "tribol::getJacobianSparseMatrix(): " << 
                 "invalid CouplingScheme id.");

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
         SLIC_WARNING("tribol::getJacobianSparseMatrix(): interface method does not return matrix data.");
         return 1;
      }
   }

} // end getJacobianSparseMatrix()

//------------------------------------------------------------------------------
int getJacobianCSRMatrix( int** I, int** J, RealT** vals, IndexT cs_id,
                          int* n_offsets, int* n_nonzero )
{
   // check to make sure input pointers are null
   if ( *I != nullptr || *J != nullptr || *vals != nullptr )
   {
      SLIC_WARNING("tribol::getJacobianCSRMatrix(): input pointers must be null.");
      return 1;
   }

   // get access to coupling scheme
   auto couplingScheme = CouplingSchemeManager::getInstance().findData(cs_id);

   // Note, SLIC_<>_ROOT macros are not here because it's possible not all ranks will have 
   // method data.
   SLIC_ERROR_IF(!couplingScheme, "tribol::getJacobianCSRMatrix(): invalid " << 
                 "CouplingScheme id.");

   switch (couplingScheme->getContactMethod())
   {
      case ALIGNED_MORTAR:
      {
         SLIC_WARNING("tribol::getJacobianCSRMatrix(): CSR format not currently implemented with " <<
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
         SLIC_WARNING("tribol::getJacobianCSRMatrix(): CSR format not currently implemented with "
                      "SINGLE_MORTAR. Use MFEM sparse matrix registration.");
         return 1;
      }
      default:
      {
         SLIC_WARNING("tribol::getJacobianCSRMatrix(): method does not return matrix data; " <<
                       "invalid call.");
         return 1;
      }
   }

} // end getJacobianCSRMatrix()

//------------------------------------------------------------------------------
int getElementBlockJacobians( IndexT cs_id, 
                              BlockSpace row_block, 
                              BlockSpace col_block,
                              const axom::Array<int>** row_elem_idx,
                              const axom::Array<int>** col_elem_idx,
                              const axom::Array<mfem::DenseMatrix>** jacobians )
{
   // get access to coupling scheme
   auto couplingScheme = CouplingSchemeManager::getInstance().findData(cs_id);

   // Note, SLIC_<>_ROOT macros are not here because it's possible not all ranks will have 
   // method data.
   SLIC_ERROR_IF(!couplingScheme, "tribol::getElementBlockJacobians(): invalid " << 
                 "CouplingScheme id.");

   SparseMode sparse_mode = couplingScheme->getEnforcementOptions().lm_implicit_options.sparse_mode;
   if (sparse_mode != SparseMode::MFEM_ELEMENT_DENSE)
   {
      SLIC_WARNING("Jacobian is assembled and can be accessed by " 
         "getMfemSparseMatrix() or getCSRMatrix(). For (unassembled) element "
         "Jacobian contributions, call setLagrangeMultiplierOptions() with "
         "SparseMode::MFEM_ELEMENT_DENSE before calling update().");
      return 1;
   }
   MethodData* method_data = couplingScheme->getMethodData();
   *row_elem_idx = &method_data->getBlockJElementIds()[static_cast<int>(row_block)];
   *col_elem_idx = &method_data->getBlockJElementIds()[static_cast<int>(col_block)];
   *jacobians = &method_data->getBlockJ()(
      static_cast<int>(row_block),
      static_cast<int>(col_block)
   );
   return 0;

} // end getElementBlockJacobians()

//------------------------------------------------------------------------------
void registerMortarGaps( IndexT mesh_id,
                         RealT * gaps )
{
   auto mesh = MeshManager::getInstance().findData(mesh_id);

   SLIC_ERROR_ROOT_IF(!mesh, "tribol::registerMortarGaps(): " << 
                      "no mesh with id " << mesh_id << " exists.");

   if (gaps == nullptr && mesh->numberOfElements() > 0)
   {
      SLIC_WARNING( "tribol::registerMortarGaps(): null pointer to gap data " << 
                    "on non-null mesh " << mesh_id << ".");
      mesh->isMeshValid() = false;
   }
   else
   {
      mesh->getNodalFields().m_node_gap = ArrayViewT<RealT>(
         gaps, mesh->numberOfNodes());
      mesh->getNodalFields().m_is_node_gap_set = true;
   }
   
} // end registerMortarGaps()

//------------------------------------------------------------------------------
void registerMortarPressures( IndexT mesh_id,
                              const RealT * pressures )
{
   auto mesh = MeshManager::getInstance().findData(mesh_id);

   SLIC_ERROR_ROOT_IF(!mesh, "tribol::registerMortarPressures(): " << 
                      "no mesh with id " << mesh_id << " exists.");

   if (pressures == nullptr && mesh->numberOfElements() > 0)
   {
      SLIC_WARNING( "tribol::registerMortarPressures(): null pointer to pressure data " << 
                    "on non-null mesh " << mesh_id << ".");
      mesh->isMeshValid() = false;
   }
   else
   {
      mesh->getNodalFields().m_node_pressure = ArrayViewT<const RealT>(
         pressures, mesh->numberOfNodes());
      mesh->getNodalFields().m_is_node_pressure_set = true;
   }
   
} // end registerMortarPressures()

//------------------------------------------------------------------------------
void registerIntNodalField( IndexT mesh_id,
                            const IntNodalFields field,
                            int * TRIBOL_UNUSED_PARAM(fieldVariable) )
{
   auto mesh = MeshManager::getInstance().findData(mesh_id);

   SLIC_ERROR_ROOT_IF(!mesh, "tribol::registerIntNodalField(): " << 
                      "no mesh with id " << mesh_id << " exists.");

   switch (field)
   {
      case UNDEFINED_INT_NODAL_FIELD:
      default:
         SLIC_ERROR_ROOT("tribol::registerIntNodalField() not yet implemented.");
   } // end switch over field

} // end registerIntNodalField()

//------------------------------------------------------------------------------
void registerRealElementField( IndexT mesh_id,
                               const RealElementFields field,
                               const RealT * fieldVariable )
{
   auto mesh = MeshManager::getInstance().findData(mesh_id);

   SLIC_ERROR_IF(!mesh, "tribol::registerRealElementField(): " << 
                 "no mesh with id " << mesh_id << " exists.");

   switch (field)
   {
      case KINEMATIC_CONSTANT_STIFFNESS:
      {
         if (fieldVariable==nullptr)
         {
            if (mesh->numberOfElements()>0)
            {
               SLIC_ERROR( "tribol::registerRealElementField(): null pointer to data for " << 
                           "'KINEMATIC_CONSTANT_STIFFNESS' on mesh " << mesh_id << ".");
               mesh->getElementData().m_is_kinematic_constant_penalty_set = false;
            }
            else
            {
               mesh->getElementData().m_is_kinematic_constant_penalty_set = true;
            }
         }
         else
         {
            mesh->getElementData().m_penalty_stiffness = *fieldVariable;
            mesh->getElementData().m_is_kinematic_constant_penalty_set = true;
         }
         break;
      }
      case RATE_CONSTANT_STIFFNESS:
      {
         if (fieldVariable==nullptr)
         {
            if (mesh->numberOfElements()>0)
            {
               SLIC_ERROR( "tribol::registerRealElementField(): null pointer to data for " << 
                           "'RATE_CONSTANT_STIFFNESS' on mesh " << mesh_id << ".");
               mesh->getElementData().m_is_rate_constant_penalty_set = false;
            }
            else
            {
               mesh->getElementData().m_is_rate_constant_penalty_set = true;
            }
         }
         else
         {
            mesh->getElementData().m_rate_penalty_stiffness = *fieldVariable;
            mesh->getElementData().m_is_rate_constant_penalty_set = true;
         }
         break;
      }
      case RATE_PERCENT_STIFFNESS:
      {
         if (fieldVariable==nullptr)
         {
            if (mesh->numberOfElements()>0)
            {
               SLIC_ERROR( "tribol::registerRealElementField(): null pointer to data for " << 
                           "'RATE_PERCENT_STIFFNESS' on mesh " << mesh_id << ".");
               mesh->getElementData().m_is_rate_percent_penalty_set = false;
            }
            else
            {
               mesh->getElementData().m_is_rate_percent_penalty_set = true;
            }
         }
         else
         {
            mesh->getElementData().m_rate_percent_stiffness = *fieldVariable;
            mesh->getElementData().m_is_rate_percent_penalty_set = true;
         }
         break;
      }
      case BULK_MODULUS:
      {
         if (fieldVariable==nullptr)
         {
            if (mesh->numberOfElements()>0)
            {
               SLIC_ERROR( "tribol::registerRealElementField(): null pointer to data for " << 
                           "'BULK_MODULUS' on mesh " << mesh_id << ".");
               mesh->getElementData().m_is_kinematic_element_penalty_set = false;
            }
            else
            {
               // set boolean to true for zero element meshes (acceptable registration)
               mesh->getElementData().m_is_kinematic_element_penalty_set = true;
            }
         }
         else
         {
            mesh->getElementData().m_mat_mod = ArrayViewT<const RealT>(fieldVariable, mesh->numberOfElements());

            // Only set boolean to true if the element thickness has been registered 
            // for nonzero element meshes. This will be true if the element thickness 
            // was registered first (need both).
            if (!mesh->getElementData().m_thickness.empty())
            {
               mesh->getElementData().m_is_kinematic_element_penalty_set = true;
            }
         }

         break;
      }
      case YOUNGS_MODULUS:
      {
         if (fieldVariable==nullptr)
         {
            if (mesh->numberOfElements()>0)
            {
               SLIC_ERROR( "tribol::registerRealElementField(): null pointer to data for " << 
                           "'YOUNGS_MODULUS' on mesh " << mesh_id << ".");
               mesh->getElementData().m_is_kinematic_element_penalty_set = false;
            }
            else
            {
               // set boolean to true for zero element meshes (acceptable registration)
               mesh->getElementData().m_is_kinematic_element_penalty_set = true;
            }
         }
         else
         {
            mesh->getElementData().m_mat_mod = ArrayViewT<const RealT>(fieldVariable, mesh->numberOfElements());

            // Only set boolean to true if the element thickness has been registered 
            // for nonzero element meshes. This will be true if the element thickness
            // was registered first (need both).
            if (!mesh->getElementData().m_thickness.empty())
            {
               mesh->getElementData().m_is_kinematic_element_penalty_set = true;
            }
         }

         break;
      }
      case ELEMENT_THICKNESS:
      {
         if (fieldVariable==nullptr)
         {
            if (mesh->numberOfElements()>0)
            {
               SLIC_ERROR( "tribol::registerRealElementField(): null pointer to data for " << 
                           "'ELEMENT_THICKNESS' on mesh " << mesh_id << ".");
               mesh->getElementData().m_is_kinematic_element_penalty_set = false;
            }
            else
            {
               // set booleans to true for zero element meshes (acceptable registration)
               mesh->getElementData().m_is_kinematic_element_penalty_set = true;
               mesh->getElementData().m_is_element_thickness_set = true;
            }
         }
         else
         {
            mesh->getElementData().m_thickness = ArrayViewT<const RealT>(fieldVariable, mesh->numberOfElements());
            mesh->getElementData().m_is_element_thickness_set = true;

            // Only set boolean to true if the material modulus has been registered for 
            // nonzero element meshes. This will set to true if the material modulus was 
            // registered first (need both).
            if (!mesh->getElementData().m_mat_mod.empty())
            {
               mesh->getElementData().m_is_kinematic_element_penalty_set = true;
            }
         }

         break;
      }
      default:
      {
         SLIC_ERROR( "tribol::registerRealElementField(): the field argument " << 
                     "on mesh " << mesh_id << " is not an accepted tribol real element field." );
      }
   } // end switch over field

} // end registerRealElementField()

//------------------------------------------------------------------------------
void registerIntElementField( IndexT mesh_id,
                              const IntElementFields field,
                              int * TRIBOL_UNUSED_PARAM(fieldVariable) )
{
   auto mesh = MeshManager::getInstance().findData(mesh_id);

   SLIC_ERROR_ROOT_IF(!mesh, "tribol::registerIntElementField(): " << 
                      "no mesh with id " << mesh_id << " exists.");

   switch (field)
   {
      case UNDEFINED_INT_ELEMENT_FIELD:
      default:
         SLIC_ERROR_ROOT("tribol::registerIntElementField() not yet implemented.");
   } // end switch over field

} // end registerIntElementField()

//------------------------------------------------------------------------------
void registerCouplingScheme( IndexT cs_id,
                             IndexT mesh_id1,
                             IndexT mesh_id2,
                             int contact_mode,
                             int contact_case,
                             int contact_method,
                             int contact_model,
                             int enforcement_method,
                             int binning_method )
{
   CouplingScheme scheme (cs_id,
                          mesh_id1,
                          mesh_id2,
                          contact_mode,
                          contact_case,
                          contact_method,
                          contact_model,
                          enforcement_method,
                          binning_method);

   // add coupling scheme to manager. Validity checks are performed in 
   // tribol::update() when each coupling scheme is initialized.
   CouplingSchemeManager::getInstance().addData(cs_id, std::move(scheme));

} // end registerCouplingScheme()

//------------------------------------------------------------------------------
void setInterfacePairs( IndexT cs_id,
                        IndexT numPairs,
                        IndexT const * const mesh_id1,
                        IndexT const * const pairType1,
                        IndexT const * const pairIndex1,
                        IndexT const * const mesh_id2,
                        IndexT const * const pairType2,
                        IndexT const * const pairIndex2 )
{
   // get access to coupling scheme
   auto couplingScheme = CouplingSchemeManager::getInstance().findData(cs_id);

   SLIC_ERROR_ROOT_IF(!couplingScheme, 
                      "tribol::setInterfacePairs(): invalid coupling scheme index.");

   auto* pairs = couplingScheme->getInterfacePairs();

   pairs->clear();

   // copy the interaction pairs
   for(int i=0; i< numPairs; ++i)
   {
      InterfacePair pair { mesh_id1[i], pairType1[i], pairIndex1[i],
                           mesh_id2[i], pairType2[i], pairIndex2[i], i };
      ContactMode mode = couplingScheme->getContactMode();

      // perform initial face-pair validity checks to add valid face-pairs 
      // to interface pair manager. Note, further computational geometry 
      // filtering will be performed on each face-pair indendifying 
      // contact candidates.
      bool check = geomFilter( pair, mode );

      if (check)
      {
         pair.isContactCandidate = true;
         pairs->addInterfacePair( pair );
      }
      else
      {
         pair.isContactCandidate = false;
      }
   }

   // Disable per-cycle rebinning
   couplingScheme->setFixedBinning(true);

} // end setInterfacePairs()

//------------------------------------------------------------------------------
int update( int cycle, RealT t, RealT &dt )
{
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
   for (auto& cs_pair : CouplingSchemeManager::getInstance())
   {
      auto& couplingScheme = cs_pair.second;

      // initialize and check for valid coupling scheme. If not valid, the coupling 
      // scheme will not be valid across all ranks and we will skip this coupling scheme
      if (!couplingScheme.init())
      {
         SLIC_WARNING_ROOT("tribol::update(): skipping invalid CouplingScheme " << 
                           cs_pair.first << "Please see warnings.");
         continue;
      }

      // perform binning between meshes on the coupling scheme
      // Note, this routine is guarded against null meshes
      couplingScheme.performBinning();

      // apply the coupling scheme. Note, there are appropriate guards against zero 
      // element meshes, or null-mesh coupling schemes
      err_cs = couplingScheme.apply( cycle, t, dt );

      if ( err_cs != 0 )
      {
         SLIC_WARNING("tribol::update(): coupling scheme " << cs_pair.first <<
                      " returned with an error.");
      }

   } // end of coupling scheme loop

   return err_cs;

} // end update()

//------------------------------------------------------------------------------
void finalize()
{
   CouplingSchemeManager::getInstance().clear();
}

//------------------------------------------------------------------------------
} // end tribol namespace

