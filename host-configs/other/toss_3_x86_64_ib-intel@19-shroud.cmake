set(CONFIG_NAME "quartz-toss_3_x86_64_ib-intel@19.0.4-fortran" CACHE PATH "") 

set(INTEL_HOME "/usr/tce/packages/intel/intel-19.0.4" CACHE PATH "")
set(CMAKE_C_COMPILER       "${INTEL_HOME}/bin/icc" CACHE PATH "")
set(CMAKE_CXX_COMPILER     "${INTEL_HOME}/bin/icpc" CACHE PATH "")
set(CMAKE_Fortran_COMPILER "${INTEL_HOME}/bin/ifort" CACHE PATH "")
set(ENABLE_FORTRAN ON CACHE BOOL "" FORCE)

set(ENABLE_MPI ON CACHE BOOL "" FORCE)
set(MPI_HOME             "/usr/tce/packages/mvapich2/mvapich2-2.3-intel-19.0.4" CACHE PATH "")
set(MPI_C_COMPILER       "${MPI_HOME}/bin/mpicc"   CACHE PATH "")
set(MPI_CXX_COMPILER     "${MPI_HOME}/bin/mpicxx"  CACHE PATH "")
set(MPI_Fortran_COMPILER "${MPI_HOME}/bin/mpifort" CACHE PATH "")
set(MPIEXEC              "/usr/bin/srun" CACHE PATH "")
set(MPIEXEC_NUMPROC_FLAG "-n" CACHE PATH "")

# Adds shroud for generating C and Fortran interfaces with "make generate" target
set(SHROUD_EXECUTABLE "/usr/apps/shroud/bin/shroud-develop" CACHE PATH "")
