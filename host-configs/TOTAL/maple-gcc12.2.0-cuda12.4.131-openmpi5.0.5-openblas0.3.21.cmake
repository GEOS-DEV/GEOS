#########################################################################
# Maple - gcc - cuda - openmpi - openblas
#########################################################################
#
# Load modules this way :
#  module purge
#  module load genesis
#  module load common
#  module load cmake/3.29.6
#  module load gcc/12.2.0
#  module load cuda/12.4.131
#  module load openmpi-gcc/5.0.5/cuda.12.4
#
#########################################################################

include(${CMAKE_CURRENT_LIST_DIR}/maple-gcc12.2.0-openmpi5.0.5-openblas0.3.21.cmake)

set( CONFIG_NAME "maple-gcc12.2.0-cuda12.4.131-openmpi5.0.5-openblas0.3.21" CACHE PATH "" )

#########################################################################
# CUDA SETUP
#########################################################################

# use :
# - CUDA 12.4.131
if( NOT DEFINED ENV{CUDA_HOME} )
    message( FATAL_ERROR "CUDA is not loaded. Please load the cuda/12.4.131 module." )
endif()

set_cuda_options()

set( ENABLE_CUDA ON CACHE BOOL "" )

#########################################################################
# OTHER OPTIONS                                                                                                                                                            
#########################################################################
# Set HYPRE to use CUDA
set( ENABLE_HYPRE_DEVICE "CUDA" CACHE PATH "" FORCE )

add_third_party_libraries()

