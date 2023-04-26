// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Tribol Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (MIT)

#ifndef SIMPLE_TRIBOL_HPP_
#define SIMPLE_TRIBOL_HPP_

#include "tribol/types.hpp"
#include "tribol/common/Parameters.hpp"

#include <string>



//------------------------------------------------------------------------------
// free functions for simple API usage
//------------------------------------------------------------------------------
/// \name Contact Library simple API methods
/// @{

/*!
 * \brief Initializes tribol and optionally initializes slic (logging library)
 *
 * \param [in] dim dimension of problem
 * \param [in] init_slic indicates if we should initialize slic
 * \return 0 if no error has occurred 
 *
 */
int Initialize(const int dim, bool init_slic = true);

/*!
 * \brief Finalizes tribol and optionally finalizes slic (logging library)
 *
 * \param [in] finalize_slic indicates if we should finalize slic
 * \return 0 if no error has occurred 
 *
 */
int Finalize(bool finalize_slic = true);


/*!
 * \brief Simple coupling setup
 *
 * \param [in] dim dimension of problem
 * \param [in] contact_method method name
 * \param [in] master_numCells the number of master surface cells
 * \param [in] master_lengthNodalData the length of the master nodal data arrays
 * \param [in] master_connectivity connectivity array for master side
 * \param [in] master_x x-coordinates of master nodes
 * \param [in] master_y y-coordinates of master nodes
 * \param [in] master_z z-coordinates of master nodes
 * \param [in] slave_numCells the number of slave surface cells
 * \param [in] slave_lengthNodalData the length of the slave nodal data arrays
 * \param [in] slave_connectivity connectivity array for slave side
 * \param [in] slave_x x-coordinates of slave nodes
 * \param [in] slave_y y-coordinates of slave nodes
 * \param [in] slave_z z-coordinates of slave nodes
 * \param [in] area_frac (optional) area fraction for overlap inclusion. Default value: 1.e-3
 * \param [in] mortar_gaps (optional) pointer to nodal mortar gap scalar field
 * \param [in] mortar_pressures (optional) pointer to nodal mortar pressures scalar field
 */
void SimpleCouplingSetup( const int dim, 
                          int contact_method,           
                          int master_numCells,
                          int master_lengthNodalData,
                          const int* master_connectivity,
                          const double* master_x,
                          const double* master_y,
                          const double* master_z,
                          int slave_numCells,
                          int slave_lengthNodalData,
                          const int* slave_connectivity,
                          const double* slave_x,
                          const double* slave_y,
                          const double* slave_z,
                          const double area_frac = 1.e-3,
                          double* mortar_gaps = nullptr,
                          double* mortar_pressures = nullptr
                        ); 

/*!
 * \brief Update per registered contact method
 *
 * \params [in,out] dt timestep
 *
 * \return 0 if no error has occurred and any update to timestep
 *
 */
int Update( double &dt);

/*!
 * \brief Gets the CSR data from the coupling scheme
 *
 * \param [out] I offsets for nonzero values array
 * \param [out] J column index array
 * \param [out] vals array of nonzero values
 * \param [out] n_offsets size of the I array
 * \param [out] n_nonzero size of the J and vals arrays
 *
 * \return 0 for success, 1 for failure
 *
 */
int GetSimpleCouplingCSR( int** I, int** J, double** vals,
                          int* n_offsets, int* n_nonzeros );

/// @}

#endif /* SIMPLE_TRIBOL_HPP_ */
