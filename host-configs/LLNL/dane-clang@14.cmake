# cmake executable path: /usr/tce/bin/cmake
cmake_minimum_required(VERSION 3.23)

set( CONFIG "clang@14")
set( CONFIG_NAME "dane-${CONFIG}" )

# Convert CMAKE_BUILD_TYPE to lowercase
string( TOLOWER "${CMAKE_BUILD_TYPE}" LOWER_CMAKE_BUILD_TYPE )

# Specify C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Set compilers
set(CMAKE_C_COMPILER /usr/tce/packages/clang/clang-14.0.6-magic/bin/clang)
set(CMAKE_CXX_COMPILER /usr/tce/packages/clang/clang-14.0.6-magic/bin/clang++)
set(MPI_C_COMPILER /usr/tce/packages/mvapich2/mvapich2-2.3.7-clang-14.0.6-magic/bin/mpicc)
set(MPI_CXX_COMPILER /usr/tce/packages/mvapich2/mvapich2-2.3.7-clang-14.0.6-magic/bin/mpicxx)

# Set the project name
project(geosx)

set( ENABLE_DOXYGEN OFF CACHE BOOL "" )
set( ENABLE_FMT OFF CACHE BOOL "" )

# Check if GEOS_PARENT_DIR is not defined by the bash script geos builder
if(NOT DEFINED GEOS_ROOT_DIR)
    # Define GEOS_ROOT_DIR with respect to this directory
    set(GEOS_ROOT_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../" CACHE PATH "")
endif()

set( GEOS_TPL_ROOT_DIR ${GEOS_ROOT_DIR}/GEOS CACHE PATH "" FORCE )
set(GEOS_TPL_DIR "${GEOS_TPL_ROOT_DIR}/tplInstall-${CONFIG_NAME}-${LOWER_CMAKE_BUILD_TYPE}" CACHE PATH "" )

# message(FATAL_ERROR "GEOS_TPL_DIR -- ${GEOS_TPL_DIR}" )

# tplInstall-dane-clang@14-release/conduit/lib/cmake/conduit/
set( CONDUIT_DIR ${GEOS_TPL_DIR}/conduit/ CACHE PATH "" FORCE )
set( PUGIXML_DIR ${GEOS_TPL_DIR}/pugixml/ CACHE PATH "" FORCE )
set( RAJA_DIR ${GEOS_TPL_DIR}/raja/ CACHE PATH "" FORCE )
set( CHAI_DIR ${GEOS_TPL_DIR}/chai/ CACHE PATH "" FORCE )
set( umpire_DIR ${GEOS_TPL_DIR}/chai/lib/cmake/umpire/ CACHE PATH "" FORCE )
set( ADIAK_DIR ${GEOS_TPL_DIR}/adiak/ CACHE PATH "" FORCE )
set( MATHPRESSO_DIR ${GEOS_TPL_DIR}/mathpresso/ CACHE PATH "" FORCE )
set( METIS_DIR ${GEOS_TPL_DIR}/metis/ CACHE PATH "" FORCE )
set( PARMETIS_DIR ${GEOS_TPL_DIR}/parmetis/ CACHE PATH "" FORCE )
set( SCOTCH_DIR ${GEOS_TPL_DIR}/scotch/ CACHE PATH "" FORCE )
set( SUPERLU_DIST_DIR ${GEOS_TPL_DIR}/superlu_dist/ CACHE PATH "" FORCE )
set( SUITESPARSE_DIR ${GEOS_TPL_DIR}/suitesparse/ CACHE PATH "" FORCE )
set( HYPRE_DIR ${GEOS_TPL_DIR}/hypre/ CACHE PATH "" FORCE )
#set( TRILINOS_DIR ${GEOS_TPL_DIR}/trilinos/lib64/cmake/Trilinos/ CACHE PATH "" FORCE )
set( TRILINOS_DIR ${GEOS_TPL_DIR}/trilinos/ CACHE PATH "" FORCE )
set( VTK_DIR ${GEOS_TPL_DIR}/vtk/ CACHE PATH "" FORCE )
set( FMT_DIR ${GEOS_TPL_DIR}/fmt/ CACHE PATH "" FORCE )
set( umpire_DIR ${GEOS_TPL_DIR}/chai/lib64/cmake/umpire/ CACHE PATH "" FORCE)

# Add a raja-config.cmake to CMAKE_PREFIX_PATH
set(CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH};${RAJA_DIR}/raja/lib/cmake/raja/")
# Add a umpire-config.cmake to CMAKE_PREFIX_PATH
set(CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH};${CHAI_DIR}/lib/cmake/umpire/")
#set(CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH};${TRILINOS_DIR}/lib64/cmake/Trilinos/")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I${PUGIXML_DIR}/include")

include(${CMAKE_CURRENT_LIST_DIR}/../../src/coreComponents/LvArray/host-configs/LLNL/dane-clang@14.cmake)

# C++
# The "-march=native -mtune=native" which LvArray adds breaks the PVT package.
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG" CACHE STRING "" FORCE)

# Fortran
set(CMAKE_Fortran_COMPILER /usr/tce/packages/gcc/gcc-12.1.1-magic/bin/gfortran CACHE PATH "")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -DNDEBUG -march=native -mtune=native" CACHE STRING "")

# MPI
set(ENABLE_MPI ON CACHE BOOL "" FORCE)
set(MPI_HOME /usr/tce/packages/mvapich2/mvapich2-2.3.7-clang-14.0.6-magic CACHE PATH "")

# Silo and VTK
set(ENABLE_SILO OFF CACHE BOOL "" FORCE)

## Set some MPM compiler options 
# set(MPMPROFILER_OPTION 1) 	# 0 off, 1 something, 2 something else.  doesn't work yet
# set(MPMP2G_DATA_OPTION 0)  	# 0 default, 1 dump P2G data size into log

#configure_file(${CMAKE_CURRENT_LIST_DIR}/../../src/coreComponents/physicsSolvers/solidMechanics/mpmConfigFiles/config.h.in ${CMAKE_CURRENT_LIST_DIR}/../../src/coreComponents/physicsSolvers/solidMechanics/mpmConfigFiles/config.h @ONLY)
#include(${CMAKE_CURRENT_LIST_DIR}/dane-base.cmake)
