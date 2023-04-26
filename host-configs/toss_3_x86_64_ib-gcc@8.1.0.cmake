set(CONFIG_NAME "quartz-toss_3_x86_64_ib-gcc@8.1.0" CACHE PATH "") 

set(CMAKE_C_COMPILER "/usr/tce/packages/gcc/gcc-8.1.0/bin/gcc" CACHE PATH "")
set(CMAKE_CXX_COMPILER "/usr/tce/packages/gcc/gcc-8.1.0/bin/g++" CACHE PATH "")

set(ENABLE_FORTRAN OFF CACHE BOOL "" FORCE)
set(ENABLE_MPI ON CACHE BOOL "" FORCE)

set(MPI_HOME             "/usr/tce/packages/mvapich2/mvapich2-2.3-gcc-8.1.0" CACHE PATH "")
set(MPI_C_COMPILER       "${MPI_HOME}/bin/mpicc"   CACHE PATH "")
set(MPI_CXX_COMPILER     "${MPI_HOME}/bin/mpicxx"  CACHE PATH "")
set(MPI_Fortran_COMPILER "${MPI_HOME}/bin/mpifort" CACHE PATH "")

set(MPIEXEC              "/usr/bin/srun" CACHE PATH "")
set(MPIEXEC_NUMPROC_FLAG "-n" CACHE PATH "")


# Temporarily disable some warnings for this configuration
set(BLT_CXX_FLAGS "-Wno-unused-parameter -Wno-unused-variable -Wno-switch -Wno-sign-compare -Wno-unused-but-set-variable -Wno-reorder" CACHE PATH "")

