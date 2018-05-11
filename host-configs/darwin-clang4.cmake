
site_name(HOST_NAME)
set(CONFIG_NAME "${HOST_NAME}-darwin-x86_64-clang@mp4" CACHE PATH "") 
message( "CONFIG_NAME = ${CONFIG_NAME}" )


set( GEOSX_TPL_ROOT_DIR "../../thirdPartyLibs/" CACHE PATH "" )


set(CMAKE_C_COMPILER "/opt/local/bin/clang-mp-4.0" CACHE PATH "")
set(CMAKE_CXX_COMPILER "/opt/local/bin/clang++-mp-4.0" CACHE PATH "")
set(ENABLE_FORTRAN OFF CACHE BOOL "" FORCE)
set(ENABLE_MPI ON CACHE PATH "")
set(MPI_C_COMPILER "/opt/local/bin/mpicc-openmpi-clang40" CACHE PATH "")
set(MPI_CXX_COMPILER "/opt/local/bin/mpicxx-openmpi-clang40" CACHE PATH "")
set(MPI_Fortran_COMPILER "/opt/local/bin/mpifort-openmpi-clang40" CACHE PATH "")
set(MPIEXEC "mpirun-openmpi-clang40" CACHE PATH "")

set(ENABLE_MATHPRESSO ON CACHE BOOL  "Enables mathpresso Plugin")


#######################################
# RAJA/CHAI SETUP
#######################################

set(CUDA_ENABLED      "OFF"       CACHE PATH "" FORCE)
set(ENABLE_OPENMP     "OFF"        CACHE PATH "" FORCE)
set(CHAI_BUILD_TYPE   "cpu-no-rm" CACHE PATH "" FORCE)
set(CHAI_ARGS         ""          CACHE PATH "" FORCE)
