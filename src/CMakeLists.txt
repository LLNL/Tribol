# Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
# other Tribol Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (MIT)

#------------------------------------------------------------------------------
# Set up tribol target
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# add tribol
#------------------------------------------------------------------------------
add_subdirectory(tribol)

#------------------------------------------------------------------------------
# add redecomp
#------------------------------------------------------------------------------
if ( TRIBOL_USE_MPI )
  add_subdirectory(redecomp)
endif()

#------------------------------------------------------------------------------
# add tribol's and redecomp's tests
#------------------------------------------------------------------------------
if ( TRIBOL_ENABLE_TESTS )
  add_subdirectory(tests)
endif()

#------------------------------------------------------------------------------
# add tribol's and redecomp's examples
#------------------------------------------------------------------------------
if ( TRIBOL_ENABLE_EXAMPLES )
  add_subdirectory(examples)
endif()
