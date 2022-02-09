include(${CMAKE_CURRENT_LIST_DIR}/../../src/coreComponents/LvArray/host-configs/ORNL/crusher-cce@13.0.1.cmake)

set( CONDUIT_DIR "${GEOSX_TPL_DIR}/conduit-0.7.2-mxjk64qwmbjzzil4c6vvlnntlririfwo/" CACHE PATH "" )
set( HDF5_DIR "${GEOSX_TPL_DIR}/hdf5-1.10.8-e36mwsprgw2fms52opc3p74hfiorhhsu/" CACHE PATH "" )
set( SILO_DIR "${GEOSX_TPL_DIR}/silo-4.10.2-gmtzhwj52xhnoodoaa6znyi2dsud33jz/" CACHE PATH "" )

set( BLAS_DIR "/opt/rocm-4.2.0/" CACHE PATH "" )
set( HYPRE_DIR "${GEOSX_TPL_DIR}/hypre-2.23.0-x44nw7jzp22jgo4lx52q2plwlh6sur37/" CACHE PATH "" )
set( SUPERLU_DIST_DIR "${GEOSX_TPL_DIR}/superlu-dist-7.1.1-iu6gvmifimd6yjyl4rnxsbpdvkkdy67u/" CACHE PATH "" )

set( PUGIXML_DIR "${GEOSX_TPL_DIR}/pugixml-1.11.4-rjfutqv6lun34dfw655kvqal6rzyy55v/" CACHE PATH "" )
set( FMT_DIR "${GEOSX_TPL_DIR}/fmt-8.0.1-c5kblzwocduybmgdccthpt254nmbyzve/" CACHE PATH "" )
set( SUITESPARSE_DIR "${GEOSX_TPL_DIR}/suite-sparse-5.10.1-jgaa5fs7xmdjfslwbyw5o3oeigw3sjvi/" CACHE PATH "" )

set(ENABLE_MATHPRESSO OFF CACHE BOOL "")

set(ENABLE_PAMELA ON CACHE BOOL "")
set(ENABLE_PVTPackage ON CACHE BOOL "")

set(ENABLE_PETSC OFF CACHE BOOL "" FORCE)

set(ENABLE_CALIPER OFF CACHE BOOL "")
set(ENABLE_PAPI OFF CACHE BOOL "")

set(ENABLE_ESSL OFF CACHE BOOL "")

set(ENABLE_TRILINOS OFF CACHE BOOL "")
set(ENABLE_VTK OFF CACHE BOOL "")

# ROCM options (blas/lapack)
set( ENABLE_ROCM ON CACHE BOOL "" FORCE )
set( ROCM_ROOT "${HIP_ROOT}" CACHE PATH "" )

# HYPRE options
set( ENABLE_HYPRE_ROCM ON CACHE BOOL "" )

set( GEOSX_BUILD_OBJ_LIBS OFF CACHE BOOL "" FORCE )
#set( BUILD_SHARED_LIBS ON CACHE BOOL "" FORCE )
#set( GEOSX_BUILD_SHARED_LIBS ON CACHE BOOL "" FORCE )

