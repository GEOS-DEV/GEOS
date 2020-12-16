site_name(HOST_NAME)
set(CONFIG_NAME "${HOST_NAME}-andre-gcc-8.4.0@ubuntu" CACHE PATH "")
message( "CONFIG_NAME = ${CONFIG_NAME}" )

set(CMAKE_C_COMPILER "/usr/bin/gcc-8" CACHE PATH "")
set(CMAKE_CXX_COMPILER "/usr/bin/g++-8" CACHE PATH "")
set(ENABLE_FORTRAN OFF CACHE BOOL "" FORCE)

set(ENABLE_MPI ON CACHE PATH "")
set(MPI_C_COMPILER "/usr/bin/mpicc" CACHE PATH "")
set(MPI_CXX_COMPILER "/usr/bin/mpicxx" CACHE PATH "")
set(MPI_Fortran_COMPILER "/usr/bin/mpifort" CACHE PATH "")

set(MPIEXEC "/usr/bin/mpirun" CACHE PATH "")
set(MPIEXEC_NUMPROC_FLAG "-n" CACHE STRING "")

set(ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "" FORCE )

set(ENABLE_PAMELA ON CACHE BOOL "" FORCE)
set(ENABLE_PVTPackage ON CACHE BOOL "" FORCE)
set(ENABLE_GEOSX_PTP ON CACHE BOOL "" FORCE)

set(ENABLE_OPENMP ON CACHE BOOL "")
set(CUDA_ENABLED OFF CACHE BOOL "")

set(ENABLE_CALIPER ON CACHE BOOL "")

set( BLAS_LIBRARIES /usr/lib/libblas.so CACHE PATH "" FORCE )
set( LAPACK_LIBRARIES /usr/lib/liblapack.so CACHE PATH "" FORCE )
set(ENABLE_PETSC OFF CACHE BOOL "Enables PETSc." FORCE)

set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG -march=native -mtune=native" CACHE STRING "")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native -mtune=native" CACHE STRING "")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -DNDEBUG -march=native -mtune=native" CACHE STRING "")
