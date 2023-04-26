// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Tribol Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (MIT)

#ifndef SRC_PHYSICS_COMMONPLANE_HPP_
#define SRC_PHYSICS_COMMONPLANE_HPP_

#include "tribol/types.hpp"
#include "Physics.hpp"

namespace tribol
{
/*!
 *
 * \brief computes penalty stiffness for Common Plane + Penalty 
 *
 * \param [in] K1 Bulk modulus for face 1
 * \param [in] t1 element thickness for face 1
 * \param [in] K2 Bulk modulus for face 2 
 * \param [in] t2 element thickness for face 2
 *
 * \return face-pair based, element-wise penalty stiffness per area 
 *         
 *
 * \pre Bulk modulus and element thickness arrays are registered by host code 
 *
 */
real ComputePenaltyStiffnessPerArea( const real K1,
                                     const real t1,
                                     const real K2,
                                     const real t2,
                                     const real tiny_length );

/*!
 *
 * \brief routine to apply interface physics in the direction normal to the interface
 *
 * \param [in] cs pointer to the coupling scheme
 *
 * \return 0 if no error
 *
 */
template< >
int ApplyNormal< COMMON_PLANE, PENALTY >( CouplingScheme const * cs );

} // end namespace tribol

#endif /* SRC_PHYSICS_COMMONPLANE_HPP_ */
