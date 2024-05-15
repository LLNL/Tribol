// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Tribol Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (MIT)

#include "examples_common.hpp" // for common functionality used in examples

#include "tribol/interface/tribol.hpp"
#include "tribol/utils/TestUtils.hpp"
#include "tribol/utils/Math.hpp"
#include "tribol/common/Parameters.hpp"
#include "tribol/mesh/MethodCouplingData.hpp"
#include "tribol/mesh/CouplingScheme.hpp"
#include "tribol/mesh/InterfacePairs.hpp"
#include "tribol/mesh/MeshData.hpp"
#include "tribol/physics/CommonPlane.hpp"
#include "tribol/geom/GeomUtilities.hpp"

#include "axom/slic.hpp"

// C/C++ includes
#include <cmath>    // for std::sin() and std::cos()
#include <string>   // for std::string and operators
#include <sstream>  // for std::ostringstream

// Example command line arguments for running this example. This test creates two rectangular blocks of dimensions (l x w x h). The dimensions are 
// set such that an initial block-intersection exists, triggering the contact interaction. The blocks are discretized per the "block#_res xx yy zz" input 
// arguments (e.g. block1_res 5 3 4, and block2_res 3 2 1), where "xx", "yy", and "zz" are the number of elements in the x, y and z directions. 
// Varying "xx" and "yy" will vary the number of contact overlaps that exist between the two surfaces, and therefore, the amount of contact work.
//
// srun -n1 ./common_plane_ex --block1_res 100 50 10  --block1_min 0 0 0  --block1_max 10 5 1    --block2_res 150 75 10  --block2_min 0 0 0.95  --block2_max 10 5 1.95
// srun -n1 ./common_plane_ex --block1_res 100 50 4   --block1_min 0 0 0. --block1_max 1 1 1.05  --block2_res 150 75 4   --block2_min 0 0 0.95  --block2_max 1 1 2
// srun -n1 ./common_plane_ex --block1_res 4 4 4      --block1_min 0 0 0. --block1_max 1 1 1.05  --block2_res 4 4 4      --block2_min 0 0 0.95  --block2_max 1 1 2

/*!
 * \brief Program main.
 *
 * This example runs the common plane + penalty algorithm for two opposing blocks
 *
 * \param [in] argc argument counter
 * \param [in] argv argument vector
 *
 * \return rc return code, a non-zero return code indicates an error.
 */
int main( int argc, char** argv )
{
  ////////////////////////////////
  //                            //
  // SETUP THE EXAMPLE AND MESH //
  //                            //
  ////////////////////////////////

  // initialize
#ifdef TRIBOL_USE_MPI
  MPI_Init( &argc, &argv );
#endif
  tribol::CommT problem_comm = TRIBOL_COMM_WORLD;
  initialize_logger( problem_comm );
 
  // instantiate test mesh object. Note, this mesh object is a Tribol 
  // utility for testing. In general, a physics application will have 
  // their own mesh data.
  tribol::TestMesh mesh;

  mesh.setupContactMeshTet( 1, 1, 1,
                            0.0, 0.0, 0.0,
                            1.0, 1.0, 1.0,
                            1, 1, 1,
                            0.0, 0.0, 1.0,
                            1.0, 1.0, 2.0,
                            0., 0. );
  
  mesh.testMeshToVtk("tet_test", 0, 0.0);

  axom::slic::flushStreams();
  finalize_logger();

#ifdef TRIBOL_USE_MPI
  MPI_Finalize();
#endif

  return 0;
}
