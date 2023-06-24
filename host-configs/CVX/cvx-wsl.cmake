site_name(HOST_NAME)
set(CONFIG_NAME "cvx-wsl" CACHE PATH "" FORCE)
message( "CONFIG_NAME=${CONFIG_NAME}" )

set(CMAKE_C_COMPILER "/usr/bin/gcc" CACHE PATH "" FORCE)
set(CMAKE_CXX_COMPILER "/usr/bin/g++" CACHE PATH "" FORCE)
set(CMAKE_Fortran_COMPILER "/usr/bin/gfortran" CACHE PATH "" FORCE)
set(ENABLE_FORTRAN OFF CACHE BOOL "" FORCE)

set(ENABLE_MPI ON CACHE PATH "" FORCE)
set(MPI_C_COMPILER "/usr/bin/mpicc" CACHE PATH "" FORCE)
set(MPI_CXX_COMPILER "/usr/bin/mpicxx" CACHE PATH "" FORCE)
set(MPI_Fortran_COMPILER "/usr/bin/mpifort" CACHE PATH "" FORCE)
set(MPIEXEC_EXECUTABLE "/usr/bin/mpirun" CACHE PATH "" FORCE)

include(${CMAKE_CURRENT_LIST_DIR}/cvx-base.cmake)
