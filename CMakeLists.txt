# Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
# other Tribol Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (MIT)

#------------------------------------------------------------------------------
# Setup Tribol project
#------------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.14)

cmake_policy(SET CMP0057 NEW)

project(tribol LANGUAGES C CXX)

if (ENABLE_FORTRAN)
   enable_language(Fortran)
endif()

#------------------------------------------------------------------------------
# Setup BLT
#------------------------------------------------------------------------------
if (DEFINED BLT_SOURCE_DIR)
    # Support having an external BLT outside of the repository
    if (NOT EXISTS ${BLT_SOURCE_DIR}/SetupBLT.cmake)
        message(FATAL_ERROR "Provided 'BLT_SOURCE_DIR' does not contain SetupBLT.cmake")
    endif()
else()
    # Use internal 'blt' submodule path if BLT_SOURCE_DIR not provided
    set(BLT_SOURCE_DIR "${PROJECT_SOURCE_DIR}/cmake/blt" CACHE PATH "")
    if (NOT EXISTS ${BLT_SOURCE_DIR}/SetupBLT.cmake)
        message(FATAL_ERROR
            "Cannot locate BLT. "
            "Either run the following two commands in your git repository: \n"
            "    git submodule init cmake/blt\n"
            "    git submodule update\n"
            "Or add -DBLT_SOURCE_DIR=/path/to/blt to your CMake command." )
    endif()
endif()

# Override a conflicting target name between blt and mfem's build systems
SET(BLT_CODE_CHECK_TARGET_NAME "code_check" CACHE STRING "")

if ("${PROJECT_SOURCE_DIR}" STREQUAL "${CMAKE_SOURCE_DIR}")
    # Set some default BLT options before loading BLT only if not included in
    # another project
    if (NOT BLT_CXX_STD)
        set(BLT_CXX_STD "c++14" CACHE STRING "")
    endif()

    set(ENABLE_ASTYLE      OFF CACHE BOOL "")
    set(ENABLE_CLANGFORMAT OFF CACHE BOOL "")
    set(ENABLE_UNCRUSTIFY  OFF CACHE BOOL "")
endif()

include(${BLT_SOURCE_DIR}/SetupBLT.cmake)

#------------------------------------------------------------------------------
# Fortran Configuration
#------------------------------------------------------------------------------
if(ENABLE_FORTRAN)
    # Check C/C++ compiler compatiblity with the Fortran compiler
    include(FortranCInterface)
    FortranCInterface_VERIFY()
    FortranCInterface_VERIFY(CXX)

    # Axom assumes that all Fortran files use free formatting
    set(CMAKE_Fortran_FORMAT FREE)
endif()

#------------------------------------------------------------------------------
# Load tribol macros and options
#------------------------------------------------------------------------------
include(cmake/TribolMacros.cmake)
include(cmake/Options.cmake)
include(cmake/TribolCompilerFlags.cmake)

#------------------------------------------------------------------------------
# Add TPLs and source directories
#------------------------------------------------------------------------------
include(cmake/SetupThirdParty.cmake)
add_subdirectory(src)

if ( TRIBOL_ENABLE_DOCS )
  add_subdirectory(docs)
endif()

#------------------------------------------------------------------------------
# Generate header file with configuration options
#------------------------------------------------------------------------------
include(cmake/TribolVersion.cmake)
include(cmake/TribolConfig.cmake)

