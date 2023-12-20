#------------------------------------------------------------------------------
# !!!! This is a generated file, edit at own risk !!!!
#------------------------------------------------------------------------------
# CMake executable path: /usr/tce/bin/cmake
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Compilers
#------------------------------------------------------------------------------
# Compiler Spec: clang@=14.0.6
#------------------------------------------------------------------------------
if(DEFINED ENV{SPACK_CC})

  set(CMAKE_C_COMPILER "/usr/WS2/smithdev/libs/serac/toss_4_x86_64_ib/2023_11_30_13_11_25/spack/lib/spack/env/clang/clang" CACHE PATH "")

  set(CMAKE_CXX_COMPILER "/usr/WS2/smithdev/libs/serac/toss_4_x86_64_ib/2023_11_30_13_11_25/spack/lib/spack/env/clang/clang++" CACHE PATH "")

  set(CMAKE_Fortran_COMPILER "/usr/WS2/smithdev/libs/serac/toss_4_x86_64_ib/2023_11_30_13_11_25/spack/lib/spack/env/clang/gfortran" CACHE PATH "")

else()

  set(CMAKE_C_COMPILER "/usr/tce/packages/clang/clang-14.0.6/bin/clang" CACHE PATH "")

  set(CMAKE_CXX_COMPILER "/usr/tce/packages/clang/clang-14.0.6/bin/clang++" CACHE PATH "")

  set(CMAKE_Fortran_COMPILER "/usr/tce/packages/gcc/gcc-10.3.1/bin/gfortran" CACHE PATH "")

endif()

#------------------------------------------------------------------------------
# MPI
#------------------------------------------------------------------------------

set(MPI_C_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.3.6-clang-14.0.6/bin/mpicc" CACHE PATH "")

set(MPI_CXX_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.3.6-clang-14.0.6/bin/mpicxx" CACHE PATH "")

set(MPI_Fortran_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.3.6-clang-14.0.6/bin/mpif90" CACHE PATH "")

set(MPIEXEC_EXECUTABLE "/usr/bin/srun" CACHE PATH "")

set(MPIEXEC_NUMPROC_FLAG "-n" CACHE STRING "")

set(ENABLE_MPI ON CACHE BOOL "")

#------------------------------------------------------------------------------
# Hardware
#------------------------------------------------------------------------------

set(ENABLE_OPENMP ON CACHE BOOL "")

#------------------------------------------------------------------------------
# Options
#------------------------------------------------------------------------------

set(BUILD_REDECOMP ON CACHE BOOL "")

#------------------------------------------------------------------------------
# TPLs
#------------------------------------------------------------------------------

set(TPL_ROOT "/usr/WS2/smithdev/libs/serac/toss_4_x86_64_ib/2023_11_30_13_11_25/clang-14.0.6" CACHE PATH "")

set(AXOM_DIR "${TPL_ROOT}/axom-0.8.1.0-uhb7sol6ykcsyhkf2py5xfkc226pqpio" CACHE PATH "")

set(MFEM_DIR "${TPL_ROOT}/mfem-4.6.1.1-o6ahxwgjlx6ousqakterfqoakydzcl5r" CACHE PATH "")

set(RAJA_DIR "${TPL_ROOT}/raja-2023.06.1-hukre2lcnxbdl3oxfr5fkunviesjr6hr" CACHE PATH "")

set(UMPIRE_DIR "${TPL_ROOT}/umpire-2023.06.0-jve7dnrevmb2jpzyo6xa2hqhdeezi4ba" CACHE PATH "")

#------------------------------------------------------------------------------
# Devtools
#------------------------------------------------------------------------------

# Root directory for generated developer tools

set(DEVTOOLS_ROOT "/usr/WS2/smithdev/devtools/toss_4_x86_64_ib/2023_11_03_16_58_50/._view/ronvl3f6draudlqfuywypbscahbjsugl" CACHE PATH "")

set(TRIBOL_ENABLE_DOCS ON CACHE BOOL "")

set(SPHINX_EXECUTABLE "${DEVTOOLS_ROOT}/python-3.10.12/bin/sphinx-build" CACHE PATH "")

set(DOXYGEN_EXECUTABLE "${DEVTOOLS_ROOT}/doxygen-1.9.6/bin/doxygen" CACHE PATH "")


