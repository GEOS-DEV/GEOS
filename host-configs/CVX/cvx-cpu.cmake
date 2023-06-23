site_name(HOST_NAME)
set(CONFIG_NAME "cvx-cpu" CACHE PATH "" FORCE)
message( "CONFIG_NAME=${CONFIG_NAME}" )

set(CMAKE_C_COMPILER "/util/gcc/gcc-9.3.0/bin/gcc" CACHE PATH "" FORCE)
set(CMAKE_CXX_COMPILER "/util/gcc/gcc-9.3.0/bin/g++" CACHE PATH "" FORCE)
set(CMAKE_Fortran_COMPILER "/vend/intel/parallel_studio_xe_2020_update2/compilers_and_libraries_2020.2.254/linux/bin/intel64/ifort" CACHE PATH "" FORCE)
set(ENABLE_FORTRAN OFF CACHE BOOL "" FORCE)

set(ENABLE_MPI ON CACHE PATH "" FORCE)
set(MPI_C_COMPILER "/vend/intel/parallel_studio_xe_2018_update3/compilers_and_libraries_2018.3.222/linux/mpi/intel64/bin/mpicc" CACHE PATH "" FORCE)
set(MPI_CXX_COMPILER "/vend/intel/parallel_studio_xe_2018_update3/compilers_and_libraries_2018.3.222/linux/mpi/intel64/bin/mpicxx" CACHE PATH "" FORCE)
set(MPI_Fortran_COMPILER "/vend/intel/parallel_studio_xe_2018_update3/compilers_and_libraries_2018.3.222/linux/mpi/intel64/bin/mpiifort" CACHE PATH "" FORCE)

set(MPIEXEC_EXECUTABLE "/vend/intel/parallel_studio_xe_2018_update3/compilers_and_libraries_2018.3.222/linux/mpi/intel64/bin/mpirun" CACHE PATH "" FORCE)

include(cvx-base.cmake)
