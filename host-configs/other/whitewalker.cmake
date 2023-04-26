set(CONFIG_NAME "rzgenie-toss_3_x86_64_ib-intel@19.0.0" CACHE PATH "") 

set(CMAKE_C_COMPILER "/usr/tce/packages/intel/intel-19.0.0/bin/icc" CACHE PATH "")
set(CMAKE_CXX_COMPILER "/usr/tce/packages/intel/intel-19.0.0/bin/icpc" CACHE PATH "")

set(ENABLE_FORTRAN OFF CACHE BOOL "" FORCE)
set(ENABLE_MPI ON CACHE BOOL "" FORCE)

set(MPI_HOME             "/usr/tce/packages/mvapich2/mvapich2-2.2-intel-19.0.0" CACHE PATH "")

set(MPI_C_COMPILER       "${MPI_HOME}/bin/mpicc"   CACHE PATH "")
set(MPI_CXX_COMPILER     "${MPI_HOME}/bin/mpicxx"  CACHE PATH "")
set(MPI_Fortran_COMPILER "${MPI_HOME}/bin/mpifort" CACHE PATH "")

set(MY_TPL_ROOT          "/usr/WS1/swopscha/mfem_codes" CACHE PATH "")

set(MFEM_DIR             "${MY_TPL_ROOT}/mfem"                         CACHE PATH "")
set(HYPRE_DIR            "${MY_TPL_ROOT}/hypre/hypre-2.11.2/src/hypre" CACHE PATH "")
set(METIS_DIR            "${MY_TPL_ROOT}/metis/parmetis-4.0.3"         CACHE PATH "")
set(PARMETIS_DIR         "${MY_TPL_ROOT}/metis/parmetis-4.0.3"         CACHE PATH "")
set(SUPERLUDIST_DIR      "${MY_TPL_ROOT}/superlu/SuperLU_DIST_5.1.3"   CACHE PATH "") 

set(MPIEXEC              "/usr/bin/srun" CACHE PATH "")
set(MPIEXEC_NUMPROC_FLAG "-n" CACHE PATH "")
