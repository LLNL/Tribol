# Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
# other Tribol Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (MIT)

# Get value for some options from corresponding blt options
cmake_dependent_option(TRIBOL_ENABLE_CUDA "Enables Tribol with CUDA support" ON "ENABLE_CUDA" OFF)
cmake_dependent_option(TRIBOL_ENABLE_HIP "Enables Tribol with HIP support" ON "ENABLE_HIP" OFF)
cmake_dependent_option(TRIBOL_USE_MPI "Enables MPI in Tribol" ON "ENABLE_MPI" OFF)
cmake_dependent_option(TRIBOL_ENABLE_TESTS "Enables Tribol Tests" ON "ENABLE_TESTS" OFF)
cmake_dependent_option(TRIBOL_ENABLE_EXAMPLES "Enables Tribol Examples" ON "ENABLE_EXAMPLES" OFF)
cmake_dependent_option(TRIBOL_ENABLE_DOCS "Enables Tribol Docs" ON "ENABLE_DOCS" OFF)

option(TRIBOL_USE_SINGLE_PRECISION "Use single-precision floating point" OFF)
option(TRIBOL_USE_64BIT_INDEXTYPE "Use 64-bit index type" OFF)

option(TRIBOL_ENABLE_SIDRE "Enables use of Sidre" OFF)
option(TRIBOL_ENABLE_SLIC "Enables use of Slic" ON)

