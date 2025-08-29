#########################################################################
# Maple - gcc - openmpi - openblas
#########################################################################
#
# Load modules this way :
#  module purge
#  module load genesis
#  module load common
#  module load cmake/3.29.6
#  module load gcc/12.2.0
#  module load openmpi-gcc/5.0.5/cuda.12.4
#
#########################################################################

include(${CMAKE_CURRENT_LIST_DIR}/maple-common.cmake)

set( CONFIG_NAME "maple-gcc12.2.0-openmpi5.0.5-openblas0.3.21" CACHE PATH "" )

#########################################################################
# COMPILER SETUP
#########################################################################

# Use :
#  - GCC
if( NOT ("$ENV{LMOD_FAMILY_COMPILER}" STREQUAL "gcc" AND "$ENV{LMOD_FAMILY_COMPILER_VERSION}" STREQUAL "12.2.0"))
    message( FATAL_ERROR "GCC is not loaded. Please load the gcc/12.2.0 module." )
endif()

set( CMAKE_C_COMPILER       "gcc"      CACHE PATH "" )
set( CMAKE_CXX_COMPILER     "g++"      CACHE PATH "" )
set( CMAKE_Fortran_COMPILER "gfortran" CACHE PATH "" )

set_compiler_options()

#########################################################################
# MPI SETUP
#########################################################################

# use :
# - OpenMPI library
if( NOT "$ENV{LMOD_FAMILY_MPI}" STREQUAL "openmpi-gcc/5.0.5" )
    message( FATAL_ERROR "OpenMPI is not loaded. Please load the openmpi-gcc/5.0.5/cuda.12.4 module." )
endif()

set( ENABLE_MPI ON CACHE BOOL "" )
set( ENABLE_OPENMP ON CACHE PATH "" FORCE )

set( MPI_C_COMPILER mpicc CACHE PATH "" FORCE )
set( MPI_CXX_COMPILER mpicxx CACHE PATH "" FORCE )
set( MPI_Fortran_COMPILER mpifort CACHE PATH "" FORCE )
set( MPIEXEC_EXECUTABLE mpiexec CACHE PATH "" FORCE )

#########################################################################
# OTHER OPTIONS                                                                                                                                                            
#########################################################################
add_third_party_libraries()
