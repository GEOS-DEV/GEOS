set(CCE_VERSION 15.0.0)

include(${CMAKE_CURRENT_LIST_DIR}/../../src/coreComponents/LvArray/host-configs/ORNL/crusher-cce@${CCE_VERSION}.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/crusher-base.cmake)

set( CONDUIT_DIR "${GEOSX_TPL_DIR}/conduit" CACHE PATH "" )
set( HDF5_DIR "${GEOSX_TPL_DIR}/hdf5" CACHE PATH "" )
set( ENABLE_SILO OFF CACHE BOOL "" )

set( BLAS_DIR "/opt/rocm-${ROCM_VERSION}/" CACHE PATH "" )

set( PUGIXML_DIR "${GEOSX_TPL_DIR}/pugixml" CACHE PATH "" )
set( FMT_DIR "${GEOSX_TPL_DIR}/fmt" CACHE PATH "" )
set( SUITESPARSE_DIR "${GEOSX_TPL_DIR}/suite-sparse" CACHE PATH "" )

# HYPRE options
set( ENABLE_HYPRE_DEVICE "CPU" CACHE STRING "" )
set( ENABLE_HYPRE_MIXINT FALSE CACHE STRING "" )
set( HYPRE_DIR "/gpfs/alpine/geo127/world-shared/hypre/hypre_v2.27.0-18-ga592bbd12_cce-15.0.0_rel" CACHE PATH "" ) # CPU SERIAL (WORKS)

# set( ENABLE_HYPRE_MIXINT TRUE CACHE STRING "" )
# set( HYPRE_DIR "/gpfs/alpine/geo127/world-shared/hypre/hypre_v2.27.0-200-ge907ce401_cce-15.0.0_mixint_rel/" CACHE PATH "" ) # CPU SERIAL MIXED INT (doesn't work on our side)

# set( HYPRE_DIR "/gpfs/alpine/geo127/world-shared/hypre/hypre_v2.27.0-18-ga592bbd12_cce-15.0.0_omp_rel" CACHE PATH "" ) # CPU OPENMP (had openmp runtime issue)

# set( ENABLE_HYPRE_DEVICE "HIP" CACHE STRING "" )
# set( HYPRE_DIR "/gpfs/alpine/geo127/world-shared/hypre/hypre_v2.27.0-200-ge907ce401_cce-15.0.0_rocm-5.4.0_umpire-2022.3.0_rel" CACHE PATH "" ) # patched hip version int32
# set( HYPRE_DIR "/gpfs/alpine/geo127/world-shared/hypre/hypre_v2.27.0-200-ge907ce401_cce-15.0.0_rocm-5.4.0_mixint_umpire-2022.3.0_rel" CACHE PATH "" ) # patched hip verison int32/64


# if( ${ENABLE_HYPRE_DEVICE} STREQUAL "HIP" )
  # set( HYPRE_DIR "/gpfs/alpine/csc326/world-shared/victorapm/hypre-v2.26.0-11-gb93beb946-cce14.0.1_rocm5.1_uvm_rel" CACHE PATH "" )
  # set( HYPRE_DIR "/gpfs/alpine/geo127/world-shared/hypre/hypre_v2.27.0-18-ga592bbd12_cce-15.0.0_rocm-5.4.0_uvm_rel" CACHE PATH "" )
  # set( HYPRE_DIR "/gpfs/alpine/geo127/world-shared/hypre/hypre_v2.27.0-18-ga592bbd12_cce-15.0.0_rocm-5.4.0_rel" CACHE PATH "" )
# else()
#   set( HYPRE_DIR "${GEOSX_TPL_DIR}/hypre" CACHE PATH "" )
# endif()

set( ENABLE_CALIPER ON CACHE BOOL "" FORCE )
set( CALIPER_DIR "${GEOSX_TPL_DIR}/caliper" CACHE PATH "" )