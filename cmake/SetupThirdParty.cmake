# Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
# other Tribol Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (MIT)

#------------------------------------------------------------------------------
# Set up tribol's TPLs
#------------------------------------------------------------------------------

message(STATUS "Configuring TPLs...\n"
               "----------------------")

set(TPL_DEPS)
include(CMakeFindDependencyMacro)
#------------------------------------------------------------------------------
# axom
#------------------------------------------------------------------------------
if (DEFINED AXOM_DIR)
  message(STATUS "Setting up external Axom TPL...")
  include(${PROJECT_SOURCE_DIR}/cmake/thirdparty/SetupAxom.cmake)

  list(APPEND TPL_DEPS axom)

else()
  message(FATAL_ERROR 
     "Axom is a required dependency for tribol. "
     "Please configure tribol with a path to axom via the AXOM_DIR variable.")
endif()


#------------------------------------------------------------------------------
# mfem
#------------------------------------------------------------------------------

if (DEFINED MFEM_DIR)
  message(STATUS "Setting up external MFEM TPL...")

  include(${PROJECT_SOURCE_DIR}/cmake/thirdparty/SetupMFEM.cmake)

  list(APPEND TPL_DEPS mfem)
else()
  message(FATAL_ERROR 
     "MFEM is a required dependency for tribol."
     " Please configure tribol with a path to axom via the MFEM_DIR variable.")
endif()


#------------------------------------------------------------------------------
# Umpire
#------------------------------------------------------------------------------

if (DEFINED UMPIRE_DIR)
  message(STATUS "Setting up external Umpire TPL...")

  include(${UMPIRE_DIR}/lib/cmake/umpire/umpire-targets.cmake)

  list(APPEND TPL_DEPS umpire)
  set(TRIBOL_USE_UMPIRE TRUE)
else()
  message(STATUS "Umpire support is OFF")
endif()

#------------------------------------------------------------------------------
# Shroud - Generates C/Fortran/Python bindings
#------------------------------------------------------------------------------
if(EXISTS ${SHROUD_EXECUTABLE})
    message(STATUS "Setting up shroud TPL...")
    execute_process(COMMAND ${SHROUD_EXECUTABLE}
                    --cmake ${CMAKE_CURRENT_BINARY_DIR}/SetupShroud.cmake
                    ERROR_VARIABLE SHROUD_cmake_error
                    RESULT_VARIABLE SHROUD_cmake_result
                    OUTPUT_STRIP_TRAILING_WHITESPACE )
    if(NOT "${SHROUD_cmake_result}" STREQUAL "0")
        message(FATAL_ERROR "Error code from Shroud: ${SHROUD_cmake_result}\n${SHROUD_cmake_error}")
    endif()

    include(${CMAKE_CURRENT_BINARY_DIR}/SetupShroud.cmake)
else()
    message(STATUS "Shroud support is OFF")
endif()

# export tribol-targets
foreach(dep ${TPL_DEPS})
  # If the target is EXPORTABLE, add it to the export set
  get_target_property(_is_imported ${dep} IMPORTED)
  if(NOT ${_is_imported})
      install(TARGETS              ${dep}
              EXPORT               tribol-targets
              DESTINATION          lib)
      # Namespace target to avoid conflicts
      set_target_properties(${dep} PROPERTIES EXPORT_NAME tribol::${dep})
  endif()
endforeach()

# export BLT targets
blt_export_tpl_targets(EXPORT tribol-targets NAMESPACE tribol)

message(STATUS "--------------------------\n"
               "Finished configuring TPLs")

