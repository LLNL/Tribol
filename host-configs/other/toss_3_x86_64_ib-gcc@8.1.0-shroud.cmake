set(CONFIG_NAME "toss_3_x86_64_ib-gcc@8.1.0-shroud" CACHE PATH "") 

set(GCC_HOME "/usr/tce/packages/gcc/gcc-8.1.0" CACHE PATH "")
set(CMAKE_C_COMPILER       "${GCC_HOME}/bin/gcc" CACHE PATH "")
set(CMAKE_CXX_COMPILER     "${GCC_HOME}/bin/g++" CACHE PATH "")
set(CMAKE_Fortran_COMPILER "${GCC_HOME}/bin/gfortran" CACHE PATH "")
set(ENABLE_FORTRAN ON CACHE BOOL "" FORCE)

set(ENABLE_MPI ON CACHE BOOL "" FORCE)
set(MPI_HOME             "/usr/tce/packages/mvapich2/mvapich2-2.3-gcc-8.1.0" CACHE PATH "")
set(MPI_C_COMPILER       "${MPI_HOME}/bin/mpicc"   CACHE PATH "")
set(MPI_CXX_COMPILER     "${MPI_HOME}/bin/mpicxx"  CACHE PATH "")
set(MPI_Fortran_COMPILER "${MPI_HOME}/bin/mpifort" CACHE PATH "")
set(MPIEXEC              "/usr/bin/srun" CACHE PATH "")
set(MPIEXEC_NUMPROC_FLAG "-n" CACHE PATH "")


# Temporarily disable some warnings for this configuration
set(BLT_CXX_FLAGS "-Wno-unused-parameter -Wno-unused-variable -Wno-switch -Wno-sign-compare -Wno-unused-but-set-variable" CACHE PATH "")

# Adds shroud for generating C and Fortran interfaces with "make generate" target
set(SHROUD_EXECUTABLE "/usr/apps/shroud/bin/shroud-develop" CACHE PATH "")
