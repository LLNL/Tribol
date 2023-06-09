#------------------------------------------------------------------------------
# !!!! This is a generated file, edit at own risk !!!!
#------------------------------------------------------------------------------
# CMake executable path: /usr/tce/packages/cmake/cmake-3.19.2/bin/cmake
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Compilers
#------------------------------------------------------------------------------
# Compiler Spec: gcc@=10.3.1
#------------------------------------------------------------------------------
if(DEFINED ENV{SPACK_CC})

  set(CMAKE_C_COMPILER "/usr/WS2/smithdev/libs/serac/toss_4_x86_64_ib/2023_05_15_16_48_01/spack/lib/spack/env/gcc/gcc" CACHE PATH "")

  set(CMAKE_CXX_COMPILER "/usr/WS2/smithdev/libs/serac/toss_4_x86_64_ib/2023_05_15_16_48_01/spack/lib/spack/env/gcc/g++" CACHE PATH "")

  set(CMAKE_Fortran_COMPILER "/usr/WS2/smithdev/libs/serac/toss_4_x86_64_ib/2023_05_15_16_48_01/spack/lib/spack/env/gcc/gfortran" CACHE PATH "")

else()

  set(CMAKE_C_COMPILER "/usr/tce/packages/gcc/gcc-10.3.1/bin/gcc" CACHE PATH "")

  set(CMAKE_CXX_COMPILER "/usr/tce/packages/gcc/gcc-10.3.1/bin/g++" CACHE PATH "")

  set(CMAKE_Fortran_COMPILER "/usr/tce/packages/gcc/gcc-10.3.1/bin/gfortran" CACHE PATH "")

endif()

#------------------------------------------------------------------------------
# MPI
#------------------------------------------------------------------------------

set(MPI_C_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.3.6-gcc-10.3.1/bin/mpicc" CACHE PATH "")

set(MPI_CXX_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.3.6-gcc-10.3.1/bin/mpicxx" CACHE PATH "")

set(MPI_Fortran_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.3.6-gcc-10.3.1/bin/mpif90" CACHE PATH "")

set(MPIEXEC_EXECUTABLE "/usr/bin/srun" CACHE PATH "")

set(MPIEXEC_NUMPROC_FLAG "-n" CACHE STRING "")

set(ENABLE_MPI ON CACHE BOOL "")

#------------------------------------------------------------------------------
# Hardware
#------------------------------------------------------------------------------

set(ENABLE_OPENMP ON CACHE BOOL "")

#------------------------------------------------------------------------------
# TPLs
#------------------------------------------------------------------------------

set(AXOM_DIR "/usr/WS2/smithdev/libs/serac/toss_4_x86_64_ib/2023_05_15_16_48_01/gcc-10.3.1/axom-0.7.0.4-xcp7mh3bxpiqylx2s5fyr7t5iqp2fddx" CACHE PATH "")

set(MFEM_DIR "/usr/WS2/smithdev/libs/serac/toss_4_x86_64_ib/2023_05_15_16_48_01/gcc-10.3.1/mfem-4.5.3.1-cw5pkjd2fblb77wng7w3jmmhblem7rsn" CACHE PATH "")

set(UMPIRE_DIR "/usr/WS2/smithdev/libs/serac/toss_4_x86_64_ib/2023_05_15_16_48_01/gcc-10.3.1/umpire-2022.03.1-r5swu32safz2wco572ge7zgpsdej46ik" CACHE PATH "")

#------------------------------------------------------------------------------
# Devtools
#------------------------------------------------------------------------------

# Root directory for generated developer tools

set(DEVTOOLS_ROOT "/usr/WS2/smithdev/devtools/toss_4_x86_64_ib/2023_05_15_14_17_13/._view/cgnz3nqe7bpcm5tl4ph5mrykvso7zsze" CACHE PATH "")

set(CLANGFORMAT_EXECUTABLE "/usr/workspace/smithdev/devtools/toss_4_x86_64_ib/latest/python-3.10.10/bin/clang-format" CACHE PATH "")

set(CLANGTIDY_EXECUTABLE "/usr/workspace/smithdev/devtools/toss_4_x86_64_ib/latest/python-3.10.10/bin/clang-tidy" CACHE PATH "")

set(ENABLE_DOCS ON CACHE BOOL "")

set(SPHINX_EXECUTABLE "${DEVTOOLS_ROOT}/python-3.10.10/bin/sphinx-build" CACHE PATH "")

set(CPPCHECK_EXECUTABLE "${DEVTOOLS_ROOT}/cppcheck-2.9/bin/cppcheck" CACHE PATH "")

set(DOXYGEN_EXECUTABLE "${DEVTOOLS_ROOT}/doxygen-1.9.6/bin/doxygen" CACHE PATH "")


