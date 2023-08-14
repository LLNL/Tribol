// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Tribol Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (MIT)

#ifndef MFEM_TRIBOL_HPP_
#define MFEM_TRIBOL_HPP_

#include "tribol/types.hpp"

#ifdef BUILD_REDECOMP

#include "tribol/common/Parameters.hpp"

// MFEM includes
#include "mfem.hpp"

namespace tribol
{

/**
 * @brief Define and register a coupling scheme over an MFEM mesh
 *
 * This method is designed to enable simple registration of a contact coupling
 * scheme over a parallel decomposed MFEM mesh.  Contact surfaces are defined
 * via a list of one or more boundary attributes, which is used to construct a
 * single contact surface mfem::ParSubMesh.  The parallel domains are then
 * rebalanced using the redecomp library and the redecomp-level mesh and fields
 * are automatically registered with Tribol.  Mesh re-decomposition using
 * redecomp is done with a call to updateMfemParallelDecomposition().
 *
 * If the registered mesh contains a higher-order Nodes grid function, then
 * low-order refinement (LOR) of the mesh is applied (if needed) to use Tribol's
 * linear contact methodologies.  LOR support is currently limited to methods
 * that do not require Jacobian calculations, and is still experimental.  The
 * low-order refinement factor equals the order of the Nodes grid function, but
 * can be overridden using setMfemLowOrderRefinedFactor().
 *
 * @param [in] cs_id Index to use for the coupling scheme
 * @param [in] mesh_id_1 The first ID of the contact surface mesh
 * @param [in] mesh_id_2 The second ID of the contact surface mesh
 * @param [in] mesh MFEM volume mesh
 * @param [in] current_coords Coordinates associated with mesh
 * @param [in] b_attributes_1 Boundary attributes defining the first mesh
 * @param [in] b_attributes_2 Boundary attributes defining the second mesh
 * @param [in] contact_mode 
 * @param [in] contact_case 
 * @param [in] contact_method 
 * @param [in] contact_model 
 * @param [in] enforcement_method 
 * @param [in] binning_method 
 */
void registerMfemCouplingScheme( integer cs_id,
                                 integer mesh_id_1,
                                 integer mesh_id_2,
                                 const mfem::ParMesh& mesh,
                                 const mfem::ParGridFunction& current_coords,
                                 const std::set<integer>& b_attributes_1,
                                 const std::set<integer>& b_attributes_2,
                                 ContactMode contact_mode,
                                 ContactCase contact_case,
                                 ContactMethod contact_method,
                                 ContactModel contact_model,
                                 EnforcementMethod enforcement_method,
                                 BinningMethod binning_method = DEFAULT_BINNING_METHOD );

/**
 * @brief Sets factor of refinement in low order refined (LOR) representation of
 * the contact mesh
 *
 * This method sets the low-order mesh refinement factor for coupling schemes
 * registered with MFEM meshes.  Low-order refinement (LOR) enables
 * representation of higher-order contact surfaces and field quantities on a
 * refined, low-order mesh.  A LOR refinement factor of e.g. 2 results in twice
 * the number of elements on the LOR mesh in each dimension.  By default, the
 * LOR refinement factor equals the order of the mesh, such that the number of
 * degrees of freedom remain constant in the higher-order mesh and the low order
 * mesh.  This method allows the LOR mesh to be further refined, if desired.
 *
 * Transfer of field quantities between higher-order and LOR representations is
 * mass preserving and is accomplished through L2 minimization.  Given a
 * higher-order MFEM mesh, LOR enables linear mesh contact methodologies to be
 * applied on higher-order meshes.
 *
 * @pre Coupling scheme cs_id must be registered using
 * registerMfemCouplingScheme()
 *
 * @param cs_id The ID of the coupling scheme
 * @param lor_factor The refinement factor of the mesh
 */
void setMfemLowOrderRefinedFactor( integer cs_id,
                                   integer lor_factor );

/**
 * @brief Registers a velocity field on a MFEM mesh-defined coupling scheme
 *
 * @pre Coupling scheme cs_id must be registered using
 * registerMfemCouplingScheme()
 *
 * @param [in] cs_id The ID of the coupling scheme with the MFEM mesh
 * @param [in] v MFEM velocity ParGridFunction defined over the parent mesh
 */
void registerMfemVelocity( integer cs_id, const mfem::ParGridFunction& v );

/**
 * @brief Returns the response (RHS) vector to a given mfem::Vector
 *
 * @pre Coupling scheme cs_id must be registered using
 * registerMfemCouplingScheme()
 * @pre Redecomp mesh must be created and up to date by calling
 * updateMfemParallelDecomposition()
 * @pre Tribol data must be up to date for current geometry by calling update()
 *
 * @param [in] cs_id The ID of the coupling scheme with the MFEM mesh
 * @param [out] r mfem::Vector of the response (RHS) vector (properly sized,
 * pre-allocated, and initialized)
 */
void getMfemResponse( integer cs_id, mfem::Vector& r );

/**
 * @brief Get assembled contact contributions for the Jacobian matrix
 *
 * @pre Coupling scheme cs_id must be registered using
 * registerMfemCouplingScheme()
 * @pre Redecomp mesh must be created and up to date by calling
 * updateMfemParallelDecomposition()
 * @pre Tribol data must be up to date for current geometry by calling update()
 *
 * This method requires registration of an mfem::ParMesh with the coupling
 * scheme. The Jacobian contributions are split into a 2x2 block structure. The
 * first row (or column) block is associated with the coordinate grid function
 * degrees of freedom and the second row (or column) block is associated with
 * the Lagrange multiplier degrees of freedom.
 *
 * @param csId Coupling scheme id with a registered MFEM mesh
 * @return Jacobian contributions as an mfem::BlockOperator
 */
std::unique_ptr<mfem::BlockOperator> getMfemBlockJacobian( integer csId );

/**
 * @brief Returns gap vector to a given mfem::Vector
 *
 * @pre Coupling scheme cs_id must be registered using
 * registerMfemCouplingScheme()
 * @pre Redecomp mesh must be created and up to date by calling
 * updateMfemParallelDecomposition()
 * @pre Response vector must be up to date for current geometry by calling
 * update()
 *
 * @param [in] cs_id Coupling scheme id with a registered MFEM mesh
 * @param [out] g Nodal gap values (values do not have to be pre-allocated) on
 * the parent-linked boundary submesh
 */
void getMfemGap( integer cs_id, mfem::Vector& g );

/**
 * @brief Returns reference to nodal pressure vector on the submesh surface
 *
 * @pre Coupling scheme cs_id must be registered using
 * registerMfemCouplingScheme()
 *
 * @param cs_id Coupling scheme id with a registered MFEM mesh 
 * @return mfem::ParGridFunction& Nodal pressure vector defined on the
 * parent-linked boundary submesh
 */
mfem::ParGridFunction& getMfemPressure( integer cs_id );

/**
 * @brief Updates mesh parallel decomposition and related grid
 * functions/Jacobian when coordinates are updated
 *
 * @pre Coupling schemes must be registered using registerMfemCouplingScheme()
 */
void updateMfemParallelDecomposition();

} /* namespace tribol */

#endif /* BUILD_REDECOMP */

#endif /* MFEM_TRIBOL_HPP_ */