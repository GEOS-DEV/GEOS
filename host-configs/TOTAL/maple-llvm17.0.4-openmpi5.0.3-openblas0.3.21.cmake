#########################################################################
# Maple - llvm - openmpi - openblas
#########################################################################
#
# Load modules this way :
#  module purge
#  module load genesis
#  module load common
#  module load cmake/3.29.6
#  module load llvm/17.0.4
#  module load openmpi-llvm/cuda.12.4
#
#########################################################################

include(${CMAKE_CURRENT_LIST_DIR}/maple-common.cmake)

set( CONFIG_NAME "maple-llvm17.0.4-openmpi5.0.3-openblas0.3.21" CACHE PATH "" )

#########################################################################
# COMPILER SETUP
#########################################################################

# Use :
#  - LLVM
if( NOT "$ENV{LMOD_FAMILY_COMPILER}" STREQUAL "llvm" )
    message( FATAL_ERROR "LLVM is not loaded. Please load the llvm/17.0.4 module." )
endif()

set( CMAKE_C_COMPILER       "clang"     CACHE PATH "" )
set( CMAKE_CXX_COMPILER     "clang++"   CACHE PATH "" )
set( CMAKE_Fortran_COMPILER "flang-new" CACHE PATH "" )

set_compiler_options()

#########################################################################
# MPI SETUP
#########################################################################

# use :
# - OpenMPI library
if( NOT "$ENV{LMOD_FAMILY_MPI}" STREQUAL "openmpi-llvm/cuda.12.4" )
    message( FATAL_ERROR "OpenMPI is not loaded. Please load the openmpi-llvm/cuda.12.4 module." )
endif()

set( ENABLE_MPI ON CACHE BOOL "" )

set( MPI_C_COMPILER mpicc CACHE PATH "" FORCE )
set( MPI_CXX_COMPILER mpicxx CACHE PATH "" FORCE )
set( MPI_Fortran_COMPILER mpifort CACHE PATH "" FORCE )
set( MPIEXEC_EXECUTABLE mpiexec CACHE PATH "" FORCE )

set( ENABLE_OPENMP ON CACHE PATH "" FORCE )

#########################################################################
# OTHER OPTIONS                                                                                                                                                            
#########################################################################
add_third_party_libraries()
