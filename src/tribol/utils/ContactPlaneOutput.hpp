// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Tribol Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (MIT)

#ifndef SRC_UTILS_CONTACTPLANEOUTPUT_HPP_
#define SRC_UTILS_CONTACTPLANEOUTPUT_HPP_

#include "tribol/types.hpp"

// AXOM includes
#include "axom/slic.hpp" 

// C++ includes
#include <string>

namespace tribol
{

// forward declarations
class ContactPlaneManager;

/*!
 *
 * \brief Write the interface mesh to vtk
 *
 * \param [in] dir output directory
 * \param [in] v_type visualization type
 * \param [in] csId coupling scheme id
 * \param [in] meshId1 id for first mesh in interface
 * \param [in] meshId2 id for second mesh in interface
 * \param [in] dim spatial dimension
 * \param [in] cycle simulation cycle
 * \param [in] t simulation time step
 *
 */
void WriteContactPlaneMeshToVtk( const std::string& dir, const VisType v_type, 
                                 const integer csId, const integer meshId1, 
                                 const integer meshId2, const integer dim,
                                 const integer cycle, const real t );


} // end of namespace "tribol"

#endif /* SRC_UTILS_CONTACTPLANEOUTPUT_HPP_ */
