# Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
# other Tribol Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (MIT)

## examples
set( contact_examples
     common_plane.cpp
     mortar_lm_patch_test.cpp
     )


set(example_depends axom tribol axom::cli11 axom::fmt)
blt_list_append(TO example_depends ELEMENTS cuda IF ENABLE_CUDA )

foreach( example ${contact_examples} )

  get_filename_component( example_name ${example} NAME_WE )

  blt_add_executable(
      NAME       ${example_name}_ex
      SOURCES    ${example}
      OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY}
      DEPENDS_ON ${example_depends}
      )


  if (ENABLE_CUDA)
    set_target_properties(${example_name}_ex PROPERTIES CUDA_SEPARABLE_COMPILATION On)
  endif()

endforeach()

## Run examples as tests, when testing enabled
if ( TRIBOL_ENABLE_TESTS AND ENABLE_TESTS )

  foreach( example common_plane mortar_lm_patch_test)
    if(TRIBOL_USE_MPI)
      blt_add_test( NAME    tribol_${example}_ex
                    COMMAND ${example}_ex 
                    NUM_MPI_TASKS 1)
    else()
      blt_add_test( NAME    tribol_${example}_ex
                    COMMAND ${example}_ex )
    endif()
  endforeach() 

endif()

# Define fortran examples, for now, this requires MPI
if(ENABLE_FORTRAN AND TRIBOL_USE_MPI)

  blt_add_executable(
    NAME       simple_tribol_interface_ex
    SOURCES    simple_tribol_interface.F
    OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY}
    DEPENDS_ON tribol
    )

  if (ENABLE_CUDA)
    set_target_properties(simple_tribol_interface_ex PROPERTIES CUDA_SEPARABLE_COMPILATION On)
  endif()

endif()

## Redecomp examples
if ( TRIBOL_USE_MPI )

  set( redecomp_examples
      domain_redecomp.cpp
      multidomain_redecomp.cpp
      matrix_redecomp.cpp
      )

  foreach( example ${redecomp_examples} )

    get_filename_component( example_name ${example} NAME_WE )

    blt_add_executable(
        NAME       ${example_name}_ex
        SOURCES    ${example}
        OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY}
        DEPENDS_ON axom mfem mpi redecomp tribol
        )

    if (ENABLE_CUDA)
      set_target_properties(${example_name}_ex PROPERTIES CUDA_SEPARABLE_COMPILATION On)
    endif()

  endforeach()
  
endif()

## Tribol + redecomp examples
if ( TRIBOL_USE_MPI )

  set( examples
      mfem_mortar_lm_patch.cpp
      mfem_common_plane.cpp
      )

  foreach( example ${examples} )
    
    get_filename_component( example_name ${example} NAME_WE )

    blt_add_executable(
        NAME       ${example_name}_ex
        SOURCES    ${example}
        OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY}
        DEPENDS_ON axom mfem mpi redecomp tribol
        )

    if (ENABLE_CUDA)
      set_target_properties(${example_name}_ex PROPERTIES CUDA_SEPARABLE_COMPILATION On)
    endif()

  endforeach()
  
endif()
