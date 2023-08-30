set(CCE_VERSION 15.0.0)

include(${CMAKE_CURRENT_LIST_DIR}/../../src/coreComponents/LvArray/host-configs/ORNL/frontier-hdf5-cce@${CCE_VERSION}.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/frontier-base.cmake)

set( CONDUIT_DIR "${GEOSX_TPL_DIR}/conduit-0.8.7" CACHE PATH "" )
#set( HDF5_DIR "${GEOSX_TPL_DIR}/hdf5-1.12.2" CACHE PATH "" )
set( HDF5_DIR "/opt/cray/pe/hdf5-parallel/1.12.2.3/crayclang/14.0/" CACHE PATH "" )

set(ENABLE_SILO FALSE CACHE BOOL "" )
# set(ENABLE_SILO TRUE CACHE BOOL "" )
# set(SILO_DIR "${GEOSX_TPL_DIR}/silo" CACHE PATH "" )

set(ENABLE_VTK FALSE CACHE BOOL "" )
# set(ENABLE_VTK TRUE CACHE BOOL "" )
# set(VTK_DIR "${GEOSX_TPL_DIR}/vtk" CACHE PATH "" )

set( BLAS_DIR "/opt/rocm-${ROCM_VERSION}/" CACHE PATH "" )

set( PUGIXML_DIR "${GEOSX_TPL_DIR}/pugixml-1.11.4" CACHE PATH "" )
set( FMT_DIR "${GEOSX_TPL_DIR}/fmt-8.0.1" CACHE PATH "" )
set( SUITESPARSE_DIR "${GEOSX_TPL_DIR}/suite-sparse-5.10.1" CACHE PATH "" )

# HYPRE options
# set( ENABLE_HYPRE_DEVICE "CPU" CACHE STRING "" )
set( ENABLE_HYPRE_DEVICE "HIP" CACHE STRING "" )

#set( ENABLE_HYPRE_MIXINT FALSE CACHE STRING "" )
#set( HYPRE_DIR "/lustre/orion/geo127/world-shared/hypre/install-int32-gpu/" CACHE PATH "" )

set( ENABLE_HYPRE_MIXINT TRUE CACHE STRING "" )
set( HYPRE_DIR "/lustre/orion/geo143/proj-shared/victorapm/projects/install/hypre/hypre_v2.29.0-9-g5d445ffc6_cce-15.0.0_rocm-5.4.0_mixint_rel" CACHE PATH "" )
#set( HYPRE_DIR "/lustre/orion/geo127/world-shared/hypre/install-mixint-gpu/" CACHE PATH "" )

set( ENABLE_CALIPER ON CACHE BOOL "" FORCE )
set( CALIPER_DIR "${GEOSX_TPL_DIR}/caliper-2.8.0" CACHE PATH "" )
