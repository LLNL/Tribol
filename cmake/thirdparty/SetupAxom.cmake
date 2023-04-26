# Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
# other Tribol Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (MIT)

#------------------------------------------------------------------------------
# Set up axom TPLs
#------------------------------------------------------------------------------

# Check paths
tribol_assert_path_exists( ${AXOM_DIR} )

find_package(axom REQUIRED
             NO_DEFAULT_PATH 
             PATHS ${AXOM_DIR}/lib/cmake)

#[==[
blt_print_target_properties(TARGET axom)
blt_print_target_properties(TARGET cli11)
blt_print_target_properties(TARGET fmt)
blt_print_target_properties(TARGET conduit::conduit)
#]==]

