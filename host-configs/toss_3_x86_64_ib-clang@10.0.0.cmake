set(CONFIG_NAME "quartz-toss_3_x86_64_ib-clang@10.0.0" CACHE PATH "") 

set(CMAKE_C_COMPILER "/usr/tce/packages/clang/clang-10.0.0/bin/clang" CACHE PATH "")
set(CMAKE_CXX_COMPILER "/usr/tce/packages/clang/clang-10.0.0/bin/clang++" CACHE PATH "")

set(ENABLE_FORTRAN OFF CACHE BOOL "" FORCE)
set(ENABLE_MPI ON CACHE BOOL "" FORCE)

set(MPI_HOME             "/usr/tce/packages/mvapich2/mvapich2-2.3-clang-10.0.0" CACHE PATH "")
set(MPI_C_COMPILER       "${MPI_HOME}/bin/mpicc"   CACHE PATH "")
set(MPI_CXX_COMPILER     "${MPI_HOME}/bin/mpicxx"  CACHE PATH "")

set(MPIEXEC              "/usr/bin/srun" CACHE PATH "")
set(MPIEXEC_NUMPROC_FLAG "-n" CACHE PATH "")

set(ENABLE_DOCS OFF CACHE BOOL "")

