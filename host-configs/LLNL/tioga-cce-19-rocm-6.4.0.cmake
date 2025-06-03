include(${CMAKE_CURRENT_LIST_DIR}/../../src/coreComponents/LvArray/host-configs/LLNL/tioga-cce-19-rocm-6.4.0.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/amdgpu-base.cmake)

# MPI
set(MPI_HOME /opt/cray/pe/mpich/8.1.33.1/ofi/crayclang/18.0 CACHE PATH "")

# GPU-aware MPI option
if( ENABLE_HYPRE_GPU_AWARE_MPI )
  set( CMAKE_HIP_LINK_FLAGS "${CMAKE_HIP_LINK_FLAGS} $ENV{PE_MPICH_GTL_DIR_amd_gfx90a} $ENV{PE_MPICH_GTL_LIBS_amd_gfx90a}" CACHE STRING "" FORCE )
endif()
