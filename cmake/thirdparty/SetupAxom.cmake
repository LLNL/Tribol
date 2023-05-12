# Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
# other Tribol Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (MIT)

#------------------------------------------------------------------------------
# Set up axom TPLs
#------------------------------------------------------------------------------

# Check paths
tribol_assert_path_exists( ${AXOM_DIR} )

find_depenency(axom REQUIRED
               NO_DEFAULT_PATH 
               PATHS ${AXOM_DIR}/lib/cmake)

# Backwards compatibility
if (NOT TARGET axom::cli11 AND TARGET cli11)
    add_library(axom::cli11 ALIAS cli11)
endif()
if (NOT TARGET axom::fmt AND TARGET fmt)
    add_library(axom::fmt ALIAS fmt)
endif()

#[==[
blt_print_target_properties(TARGET axom)
blt_print_target_properties(TARGET axom::cli11)
blt_print_target_properties(TARGET axom::fmt)
blt_print_target_properties(TARGET conduit::conduit)
#]==]

