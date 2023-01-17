include(${CMAKE_CURRENT_LIST_DIR}/../../src/coreComponents/LvArray/host-configs/ORNL/crusher-cce@14.0.2.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/crusher-base.cmake)

set( CONDUIT_DIR "${GEOSX_TPL_DIR}/conduit-0.7.2-do5xtin3g5k4wifbnt5gayp7mpasiyz7/" CACHE PATH "" )
set( HDF5_DIR "${GEOSX_TPL_DIR}/hdf5-1.12.1-ogpok7wpfhdmmhdfrjqcc6udqn72zdkp/" CACHE PATH "" )
set( SILO_DIR "${GEOSX_TPL_DIR}/silo-4.11-acztirtk7r2bh5ywg6mpau6ytdeklaiw/" CACHE PATH "" )

set( BLAS_DIR "/opt/rocm-5.2.0/" CACHE PATH "" )

set( PUGIXML_DIR "${GEOSX_TPL_DIR}/pugixml-1.11.4-ppqkqluwqwxeyglelnpbcznigaw5pzyh/" CACHE PATH "" )
set( FMT_DIR "${GEOSX_TPL_DIR}/fmt-8.0.1-aaxwzvbmpej3aqfxwwxd36fmen7kd3z5/" CACHE PATH "" )
set( SUITESPARSE_DIR "${GEOSX_TPL_DIR}/suite-sparse-5.10.1-v6tr6m5ykekz2dtdn2ahara42nynzzbv/" CACHE PATH "" )

# HYPRE options
set( ENABLE_HYPRE_DEVICE "HIP" CACHE BOOL "" )

if( ${ENABLE_HYPRE_DEVICE} STREQUAL "HIP" )
  set( HYPRE_DIR "/gpfs/alpine/csc326/world-shared/victorapm/hypre-b378dca-cce14.0.1_rocm5.1_uvm_rel" CACHE PATH "" ) # victor's debug build
  # set( HYPRE_DIR "${GEOSX_TPL_DIR}/hypre-develop-yyketmjppcbcqjmf2z5q5j6jhee52sog" CACHE PATH "" ) # dev debug build
  # set( HYPRE_DIR "${GEOSX_TPL_DIR}/hypre-develop-jsr2pkhycafvztrl6s5dj2t3temd6q3f" CACHE PATH "" ) # debug build
  # set( HYPRE_DIR "${GEOSX_TPL_DIR}/hypre-develop-tumd2xi7jd5fbcb5rtvlydtssdje3h2v" CACHE PATH "" ) # unified memory
  # set( HYPRE_DIR "${GEOSX_TPL_DIR}/hypre-develop-" CACHE PATH "" ) # w/out UM
else()
  # set( HYPRE_DIR "${GEOSX_TPL_DIR}/hypre-develop-26xg3xxx6q3q3ineijpty36eorhsp7mh" CACHE PATH "" )
endif()

set( SUPERLU_DIST_DIR "/sw/crusher/spack-envs/base/opt/cray-sles15-zen3/cce-14.0.0/superlu-dist-7.2.0-kwzr52u5mhumb6exihw2did6h7zkiksq" CACHE PATH "" )
