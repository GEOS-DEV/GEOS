site_name(HOST_NAME)
set(CONFIG_NAME "cvx-cpu" CACHE PATH "" FORCE)
message( "CONFIG_NAME=${CONFIG_NAME}" )

set(GCC_ROOT "/util/gcc/gcc-9.3.0" CACHE PATH "")
include(${CMAKE_CURRENT_LIST_DIR}/cvx-gcc.cmake)

set(CMAKE_Fortran_COMPILER "/vend/intel/parallel_studio_xe_2020_update2/compilers_and_libraries_2020.2.254/linux/bin/intel64/ifort" CACHE PATH "" FORCE)

set(MPI_ROOT "/vend/intel/parallel_studio_xe_2018_update3/compilers_and_libraries_2018.3.222/linux/mpi/intel64" CACHE PATH "" FORCE)
include(${CMAKE_CURRENT_LIST_DIR}/cvx-mpi.cmake)

include(${CMAKE_CURRENT_LIST_DIR}/cvx-base.cmake)
