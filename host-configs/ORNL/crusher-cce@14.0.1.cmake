include(${CMAKE_CURRENT_LIST_DIR}/../../src/coreComponents/LvArray/host-configs/ORNL/crusher-cce@14.0.1.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/crusher-base.cmake)

set( CONDUIT_DIR "${GEOSX_TPL_DIR}/conduit-0.7.2-nwavjccfppg7a52fjsyewniofvhf75vu/" CACHE PATH "" )
set( HDF5_DIR "${GEOSX_TPL_DIR}/hdf5-1.12.1-vmtu5oxldtyuwatgxeigk76k7jhrvnp6/" CACHE PATH "" )
set( SILO_DIR "${GEOSX_TPL_DIR}/silo-4.11-gafm2b2huoks57k7cni7q64qlzbwvj32/" CACHE PATH "" )

set( BLAS_DIR "/opt/rocm-5.1.0/" CACHE PATH "" )

set( PUGIXML_DIR "${GEOSX_TPL_DIR}/pugixml-1.11.4-vcc5zgr3nfifb4kb6q3pkf25vfutanui/" CACHE PATH "" )
set( FMT_DIR "${GEOSX_TPL_DIR}/fmt-8.0.1-pjcagpw5nphzst552jgir6kqlnosc35u" CACHE PATH "" )
set( SUITESPARSE_DIR "${GEOSX_TPL_DIR}/suite-sparse-5.10.1-vouqh5scdgojuqw2fkdr5idnmbuyoeqq/" CACHE PATH "" )

# HYPRE options
set( ENABLE_HYPRE_DEVICE "HIP" CACHE BOOL "" )
if( ${ENABLE_HYPRE_DEVICE} STREQUAL "HIP" )
  set( HYPRE_DIR "/gpfs/alpine/csc326/world-shared/victorapm/hypre-b378dca-cce14.0.1_rocm5.1_uvm_rel" CACHE PATH "" ) # victor's build
else()
  set( HYPRE_DIR "${GEOSX_TPL_DIR}/hypre-develop-26xg3xxx6q3q3ineijpty36eorhsp7mh" CACHE PATH "" )
endif()

set( ENABLE_CALIPER ON CACHE BOOL "" FORCE )
#set( ADIAK_DIR "${GEOSX_TPL_DIR}/adiak-0.2.1-zvvjdmdwh3x2wijmwferzyeunujmep46" CACHE PATH "" )
set( CALIPER_DIR "${GEOSX_TPL_DIR}/caliper-2.8.0-ux2oi5h55i5j7fokzujjst5uwhw2yfeo" CACHE PATH "" )