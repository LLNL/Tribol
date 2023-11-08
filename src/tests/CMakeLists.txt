# Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
# other Tribol Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (MIT)

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Add single source tests
#------------------------------------------------------------------------------
set( tribol_tests
     tribol_check_tpl.cpp
     tribol_common_plane_penalty.cpp
     tribol_common_plane_gap_rate.cpp
     tribol_comp_geom.cpp
     tribol_coupling_scheme.cpp
     tribol_coupling_scheme_manager.cpp
     tribol_enforcement_options.cpp
     tribol_hex_mesh.cpp
     tribol_inv_iso.cpp
     tribol_iso_integ.cpp
     tribol_math.cpp
     tribol_mortar_data_geom.cpp
     tribol_mortar_data_weights.cpp
     tribol_mortar_force.cpp
     tribol_mortar_gap.cpp
     tribol_mortar_jacobian.cpp
     tribol_mortar_lm_patch_test.cpp
     tribol_mortar_sparse_weights.cpp
     tribol_mortar_wts.cpp
     tribol_nodal_nrmls.cpp
     tribol_quad_integ.cpp
     tribol_tet_mesh.cpp
     tribol_timestep_vote.cpp
     tribol_twb_integ.cpp
     )

set(test_depends axom mfem tribol gtest)
blt_list_append(TO test_depends ELEMENTS cuda IF ENABLE_CUDA )

foreach( test ${tribol_tests} )

  get_filename_component( test_name ${test} NAME_WE )

  blt_add_executable(
    NAME ${test_name}_test
    SOURCES ${test}
    OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
    DEPENDS_ON ${test_depends}
    FOLDER tribol/tests )

  blt_add_test( NAME ${test_name}
                COMMAND ${test_name}_test )

  if (ENABLE_CUDA)
    set_target_properties(${test_name}_test PROPERTIES CUDA_SEPARABLE_COMPILATION On)
  endif()

endforeach()

#------------------------------------------------------------------------------
# Add redecomp tests
#------------------------------------------------------------------------------
if ( BUILD_REDECOMP AND TRIBOL_USE_MPI )

  set( redecomp_tests
      redecomp_transfer.cpp
      redecomp_multitransfer.cpp
      redecomp_massmatrix.cpp
      redecomp_rectmatrix.cpp
      redecomp_sparsematrix.cpp
      )

  set(redecomp_test_depends axom mfem redecomp tribol gtest)
  blt_list_append(TO redecomp_test_depends ELEMENTS cuda IF ENABLE_CUDA )

  foreach( test ${redecomp_tests} )

    get_filename_component( test_name ${test} NAME_WE )

    blt_add_executable(
      NAME ${test_name}_test
      SOURCES ${test}
      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
      DEPENDS_ON ${redecomp_test_depends}
      FOLDER redecomp/tests )

    foreach ( nranks RANGE 1 3 )
      blt_add_test( NAME ${test_name}_${nranks}rank
                    COMMAND ${test_name}_test
                    NUM_MPI_RANKS ${nranks} )
    endforeach()

    if (ENABLE_CUDA)
      set_target_properties(${test_name}_test PROPERTIES CUDA_SEPARABLE_COMPILATION On)
    endif()

  endforeach()

endif()

#------------------------------------------------------------------------------
# Add tribol + redecomp tests
#------------------------------------------------------------------------------
if ( BUILD_REDECOMP AND TRIBOL_USE_MPI )

  set( combined_tests
      tribol_mfem_common_plane.cpp
      tribol_mfem_mortar_lm.cpp
      )

  set(combined_test_depends axom mfem redecomp tribol gtest)
  blt_list_append(TO combined_test_depends ELEMENTS cuda IF ENABLE_CUDA )

  foreach( test ${combined_tests} )

    get_filename_component( test_name ${test} NAME_WE )

    blt_add_executable(
      NAME ${test_name}_test
      SOURCES ${test}
      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
      DEPENDS_ON ${combined_test_depends}
      FOLDER tribol/tests )

    foreach ( nranks RANGE 1 3 )
      blt_add_test( NAME ${test_name}_${nranks}rank
                    COMMAND ${test_name}_test
                    NUM_MPI_RANKS ${nranks} )
    endforeach()

    if (ENABLE_CUDA)
      set_target_properties(${test_name}_test PROPERTIES CUDA_SEPARABLE_COMPILATION On)
    endif()

  endforeach()

endif()

#------------------------------------------------------------------------------
# Add smoke tests
#------------------------------------------------------------------------------

if( MFEM_FOUND )
  set(mfem_smoke_depends mfem gtest)
  blt_list_append(TO mfem_smoke_depends ELEMENTS umpire cuda IF ENABLE_CUDA )
  blt_add_executable(
    NAME mfem_smoke_test
    SOURCES mfem_smoke.cpp
    OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
    DEPENDS_ON ${mfem_smoke_depends}
    FOLDER tribol/tests )

  blt_add_test( NAME mfem_smoke_test
                COMMAND mfem_smoke_test )

  if (ENABLE_CUDA)
    set_target_properties(mfem_smoke_test PROPERTIES CUDA_SEPARABLE_COMPILATION On)
  endif()

endif()

