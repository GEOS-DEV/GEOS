
site_name(HOST_NAME)

set(CMAKE_C_COMPILER "/usr/bin/gcc" CACHE PATH "" FORCE)
set(CMAKE_CXX_COMPILER "/usr/bin/g++" CACHE PATH "" FORCE)
set(ENABLE_FORTRAN ON CACHE BOOL "" FORCE)

set(ENABLE_MPI ON CACHE PATH "" FORCE)
set(MPI_C_COMPILER "mpicc" CACHE PATH "" FORCE)
set(MPI_CXX_COMPILER "mpic++" CACHE PATH "" FORCE)
set(OMPI_C_COMPILER "/usr/bin/gcc" CACHE PATH "" FORCE)
set(OMPI_CXX_COMPILER "/usr/bin/g++" CACHE PATH "" FORCE)
set(MPIEXEC "mpirun" CACHE PATH "" FORCE)
set( ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "" FORCE )

set(ENABLE_PYGEOSX ON CACHE BOOL "")
set(Python3_EXECUTABLE "~/GEOSX_env/utilities/python/Python-3.8.5/build/bin/python3.8" CACHE PATH "")
set(CMAKE_CXX_FLAGS "-fno-lto make")
set(CMAKE_C_FLAGS "-fno-lto make")

set(GEOSX_TPL_DIR "~/GEOSX_env/thirdPartyLibs/install-configUbuntu21-debug" CACHE PATH "" FORCE )
include(~/GEOSX_env/GEOSX/host-configs/tpls.cmake)
