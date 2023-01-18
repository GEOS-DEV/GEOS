include(${CMAKE_CURRENT_LIST_DIR}/../../src/coreComponents/LvArray/host-configs/LLNL/tioga-cce@14.0.1.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/tioga-base.cmake)

set( CONDUIT_DIR "${GEOSX_TPL_DIR}/conduit-0.8.6-5p5sbys4qujoqi5fl42pg7hwhbmbvpvb" CACHE PATH "" )
set( HDF5_DIR "${GEOSX_TPL_DIR}/hdf5-1.12.2-e37kirkodbtq4l5j24fugc57nwo4yctv" CACHE PATH "" )

set( BLAS_DIR "/opt/rocm-5.1.0/" CACHE PATH "" )

set( PUGIXML_DIR "${GEOSX_TPL_DIR}/pugixml-1.11.4-miwfq6cfetizlzytthfjslvfftvwrkoe" CACHE PATH "" )
set( FMT_DIR "${GEOSX_TPL_DIR}/fmt-8.0.1-7ftsg6xi7spadqd22f6t4k5acnfrcm5t" CACHE PATH "" )
set( SUITESPARSE_DIR "${GEOSX_TPL_DIR}/suite-sparse-5.10.1-7sa2amcy74brjea2vodxrafs7mzhgznz" CACHE PATH "" )

# HYPRE options
set( ENABLE_HYPRE_DEVICE "HIP" CACHE STRING "" )
set( HYPRE_DIR "${GEOSX_TPL_DIR}/hypre-develop-bsuoncb2sk3o4ocjatxqt4vkpvg6q3r2" CACHE PATH "" )

set( ENABLE_CALIPER ON CACHE BOOL "" FORCE )
#set( ADIAK_DIR "${GEOSX_TPL_DIR}/adiak-0.2.1-zvvjdmdwh3x2wijmwferzyeunujmep46" CACHE PATH "" )
set( CALIPER_DIR "${GEOSX_TPL_DIR}/caliper-2.8.0-cjnthha7sceu73zoxuu73f6mu64zbmui" CACHE PATH "" )
