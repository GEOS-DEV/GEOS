include(${CMAKE_CURRENT_LIST_DIR}/../../src/coreComponents/LvArray/host-configs/LLNL/tioga-cce-15.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/tioga-base.cmake)

set( CONDUIT_DIR "${GEOSX_TPL_DIR}/conduit" CACHE PATH "" )
set( HDF5_DIR "${GEOSX_TPL_DIR}/hdf5" CACHE PATH "" )

set( BLAS_DIR "/opt/rocm-5.4.3/" CACHE PATH "" )

set( PUGIXML_DIR "${GEOSX_TPL_DIR}/pugixml" CACHE PATH "" )
set( FMT_DIR "${GEOSX_TPL_DIR}/fmt" CACHE PATH "" )
set( SUITESPARSE_DIR "${GEOSX_TPL_DIR}/suite-sparse" CACHE PATH "" )

# HYPRE options
set( ENABLE_HYPRE_DEVICE "HIP" CACHE STRING "" )
set( ENABLE_HYPRE_MIXINT FALSE CACHE STRING "" )
set( HYPRE_DIR "${GEOSX_TPL_DIR}/hypre" CACHE PATH "" )

set( ENABLE_CALIPER ON CACHE BOOL "" FORCE )
set( CALIPER_DIR "${GEOSX_TPL_DIR}/caliper" CACHE PATH "" )

# haven't build I/O TPLs on tioga yet
set( ENABLE_SILO OFF CACHE BOOL "" FORCE )