// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Tribol Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (MIT)

#include "tribol.hpp"

// Tribol includes
#include "tribol/interface/simple_tribol.hpp"
#include "tribol/common/logger.hpp"
#include "tribol/common/Parameters.hpp"
#include "tribol/types.hpp"

// Axom includes
#include "axom/core.hpp"
#include "axom/slic.hpp"

// C/C++ includes
#include <string>
#include <unordered_map>
#include <fstream>

// Axom Aliases
namespace numerics = axom::numerics;

//------------------------------------------------------------------------------
// Interface Implementation
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
// free functions for simple API usage
//------------------------------------------------------------------------------

int Initialize(const int dim, bool init_slic)
{
   // initialize slic
   if(init_slic)
   {
      axom::slic::finalize();
      axom::slic::initialize();
      std::string format = "[<LEVEL>]: <MESSAGE> \n";
      axom::slic::setLoggingMsgLevel( axom::slic::message::Info );

      axom::slic::addStreamToAllMsgLevels(
         new axom::slic::GenericOutputStream( &std::cout,format) );
   }

   // Initialize tribol
   tribol::CommType problem_comm = TRIBOL_COMM_WORLD;
   tribol::initialize( dim, problem_comm );

   return 0;
}

int Finalize(bool finalize_slic)
{
   // finalize tribol
   tribol::finalize();

   // finalize slic
   if(finalize_slic)
   {
      axom::slic::finalize();
   }
   
   return 0;
}

void SimpleCouplingSetup( const int dim, 
                          int cell_type,
                          int contact_method,           
                          int mortar_numCells,
                          int mortar_lengthNodalData,
                          const int* mortar_connectivity,
                          const double* mortar_x,
                          const double* mortar_y,
                          const double* mortar_z,
                          int nonmortar_numCells,
                          int nonmortar_lengthNodalData,
                          const int* nonmortar_connectivity,
                          const double* nonmortar_x,
                          const double* nonmortar_y,
                          const double* nonmortar_z,
                          const double area_frac,
                          double* mortar_gaps,
                          double* mortar_pressures)
{
   if (contact_method != tribol::MORTAR_WEIGHTS)
   {
      SLIC_ERROR( "SimpleCouplingSetup: simple API only works " << 
                  "for MORTAR_WEIGHTS method." );
   }

   // register mortar mesh
   int mortarMeshId = 0;
   tribol::registerMesh( mortarMeshId, mortar_numCells, 
                         mortar_lengthNodalData,
                         mortar_connectivity, cell_type,
                         mortar_x, mortar_y, mortar_z );

   // register nonmortar mesh
   int nonmortarMeshId = 1;
   tribol::registerMesh( nonmortarMeshId, nonmortar_numCells,
                         nonmortar_lengthNodalData,
                         nonmortar_connectivity, cell_type,
                         nonmortar_x, nonmortar_y, nonmortar_z );

   // Register mortar gaps and pressures, if provided
   if( mortar_gaps != nullptr)
   {
      tribol::registerRealNodalField( nonmortarMeshId, tribol::MORTAR_GAPS, mortar_gaps);
   }
   if( mortar_pressures != nullptr)
   {
      tribol::registerRealNodalField( nonmortarMeshId, tribol::MORTAR_PRESSURES, mortar_pressures);
   }

   // set contact area fraction 
   tribol::setContactAreaFrac( area_frac );
   tribol::setPlotCycleIncrement(1);

   // note the use of NULL_ENFORCEMENT reflects that this routine is used 
   // to initially setup tests for MORTAR_WEIGHTS only!
   tribol::registerCouplingScheme( 0, mortarMeshId, nonmortarMeshId,
                                   tribol::SURFACE_TO_SURFACE,
                                   tribol::AUTO,
                                   contact_method,
                                   tribol::NULL_MODEL,
                                   tribol::NULL_ENFORCEMENT );

   // set enforcement options for MORTAR_WEIGHTS
   tribol::setLagrangeMultiplierOptions( 0, tribol::ImplicitEvalMode::MORTAR_WEIGHTS_EVAL, 
                                         tribol::SparseMode::MFEM_LINKED_LIST );

   axom::slic::flushStreams();

   return;
}

//------------------------------------------------------------------------------
int Update( double &dt )
{
   int err = tribol::update( 1, 1., dt );   

   axom::slic::flushStreams();

   return err;
}

//------------------------------------------------------------------------------
int GetSimpleCouplingCSR( int** I, int** J, double** vals,
                          int* n_offsets, int* n_nonzeros )
{
   int err = tribol::getCSRMatrix( I, J, vals, 0, n_offsets, n_nonzeros );
   return err;
}

