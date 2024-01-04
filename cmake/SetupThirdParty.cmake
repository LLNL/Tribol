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
# Create global variable to toggle between GPU targets
#------------------------------------------------------------------------------
if(TRIBOL_ENABLE_CUDA)
  set(tribol_device_depends cuda CACHE STRING "" FORCE)
endif()
if(TRIBOL_ENABLE_HIP)
  set(tribol_device_depends blt::hip CACHE STRING "" FORCE)
endif()

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
     "Please configure tribol with a path to Axom via the AXOM_DIR variable.")
endif()


#------------------------------------------------------------------------------
# mfem
#------------------------------------------------------------------------------

if (TARGET mfem)
    # Case: Tribol included in project that also creates an mfem target, no need to recreate mfem
    # Note - white238: I can't seem to get this to pass install testing due to mfem being included
    # in multiple export sets
    message(STATUS "MFEM support is ON, using existing mfem target")
    # Add it to this export set but don't prefix it with tribol::
    install(TARGETS              mfem
            EXPORT               tribol-targets
            DESTINATION          lib)
    set(MFEM_FOUND TRUE CACHE BOOL "" FORCE)
elseif (DEFINED MFEM_DIR)
  message(STATUS "Setting up external MFEM TPL...")

  include(${PROJECT_SOURCE_DIR}/cmake/thirdparty/SetupMFEM.cmake)

  list(APPEND TPL_DEPS mfem)
else()
  message(FATAL_ERROR 
     "MFEM is a required dependency for tribol. "
     "Please configure tribol with a path to MFEM via the MFEM_DIR variable.")
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

#---------------------------------------------------------------------------
# Remove non-existant INTERFACE_INCLUDE_DIRECTORIES from imported targets
# to work around CMake error
#---------------------------------------------------------------------------
set(_imported_targets
    axom
    axom::mfem
    conduit
    conduit::conduit_mpi
    conduit::conduit
    conduit_relay_mpi
    conduit_relay_mpi_io
    conduit_blueprint
    conduit_blueprint_mpi)

foreach(_target ${_imported_targets})
    if(TARGET ${_target})
        message(STATUS "Removing non-existant include directories from target[${_target}]")

        get_target_property(_dirs ${_target} INTERFACE_INCLUDE_DIRECTORIES)
        set(_existing_dirs)
        foreach(_dir ${_dirs})
            if (EXISTS "${_dir}")
                list(APPEND _existing_dirs "${_dir}")
            endif()
        endforeach()
        if (_existing_dirs)
            set_target_properties(${_target} PROPERTIES
                                  INTERFACE_INCLUDE_DIRECTORIES "${_existing_dirs}" )
        endif()
    endif()
endforeach()

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

