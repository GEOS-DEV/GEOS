include(${CMAKE_CURRENT_LIST_DIR}/../../src/coreComponents/LvArray/host-configs/ORNL/spock-cce@12.0.3.cmake)

set( CONDUIT_DIR "${GEOSX_TPL_DIR}/conduit-0.7.2-qbfiuuhnlvw5blwgxmbfe2o4kpandavv/" CACHE PATH "" )
set( HDF5_DIR "${GEOSX_TPL_DIR}/hdf5-1.10.7-7pg4zt3u32shanxvygiza5gqaqrf5v3z/" CACHE PATH "" )
set( SILO_DIR "${GEOSX_TPL_DIR}/silo-4.10.2-sefw452tklfege3ivr3yytf7xypghksp/" CACHE PATH "" )

set( BLAS_DIR "/opt/rocm-4.2.0/" CACHE PATH "" )
set( HYPRE_DIR "${GEOSX_TPL_DIR}/hypre-2.23.0-fogytqv6jjuahwack4ijl6aoi3mi6i73/" CACHE PATH "" )
set( SUPERLU_DIST_DIR "${GEOSX_TPL_DIR}/superlu-dist-7.1.1-7qflwyslkjktb7eoaqqpt2ycobbtpcag/" CACHE PATH "" )

set( PUGIXML_DIR "${GEOSX_TPL_DIR}/pugixml-1.11.4-fmyqr3objxmxdunimdyfzzglrxexhepn/" CACHE PATH "" )
set( FMT_DIR "${GEOSX_TPL_DIR}/fmt-8.0.1-2kmbox3ikrgej4oousn3ba5ufoshpudt/" CACHE PATH "" )
set( SUITESPARSE_DIR "${GEOSX_TPL_DIR}/suite-sparse-5.10.1-e6hdpuyd7r4j4bpvmkots4q6nelly5hk/" CACHE PATH "" )

# asmjit doesn't work on PowerPC
set(ENABLE_MATHPRESSO OFF CACHE BOOL "")

set(ENABLE_PAMELA ON CACHE BOOL "")
set(ENABLE_PVTPackage ON CACHE BOOL "")

set(ENABLE_PETSC OFF CACHE BOOL "" FORCE)

set(ENABLE_CALIPER OFF CACHE BOOL "")
set(ENABLE_PAPI OFF CACHE BOOL "")

set(ENABLE_ESSL OFF CACHE BOOL "")

set(ENABLE_METIS OFF CACHE BOOL "")
set(ENABLE_PARMETIS OFF CACHE BOOL "")
set(ENABLE_SUPERLU_DIST OFF CACHE BOOL "")

set(ENABLE_TRILINOS OFF CACHE BOOL "")
set(ENABLE_VTK OFF CACHE BOOL "")

# ROCM options (blas/lapack)
set( ENABLE_ROCM ON CACHE BOOL "" FORCE )
set( ROCM_ROOT "${HIP_ROOT}" CACHE PATH "" )

set( ENABLE_HYPRE_ROCM ON CACHE BOOL "" )