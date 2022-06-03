include(${CMAKE_CURRENT_LIST_DIR}/../../src/coreComponents/LvArray/host-configs/ORNL/crusher-cce@14.0.0.cmake)

set( CONDUIT_DIR "${GEOSX_TPL_DIR}/conduit-0.7.2-nwavjccfppg7a52fjsyewniofvhf75vu/" CACHE PATH "" )
set( HDF5_DIR "${GEOSX_TPL_DIR}/hdf5-1.12.1-vmtu5oxldtyuwatgxeigk76k7jhrvnp6/" CACHE PATH "" )
set( SILO_DIR "${GEOSX_TPL_DIR}/silo-4.11-gafm2b2huoks57k7cni7q64qlzbwvj32/" CACHE PATH "" )

set( BLAS_DIR "/opt/rocm-5.1.0/" CACHE PATH "" )

set( PUGIXML_DIR "${GEOSX_TPL_DIR}/pugixml-1.11.4-vcc5zgr3nfifb4kb6q3pkf25vfutanui/" CACHE PATH "" )
set( FMT_DIR "${GEOSX_TPL_DIR}/fmt-8.0.1-pjcagpw5nphzst552jgir6kqlnosc35u" CACHE PATH "" )
set( SUITESPARSE_DIR "${GEOSX_TPL_DIR}/suite-sparse-5.10.1-vouqh5scdgojuqw2fkdr5idnmbuyoeqq/" CACHE PATH "" )

set( ENABLE_MATHPRESSO OFF CACHE BOOL "" )

set( ENABLE_PAMELA ON CACHE BOOL "" )
set( ENABLE_PVTPackage ON CACHE BOOL "" )

set( ENABLE_PETSC OFF CACHE BOOL "" FORCE )

set( ENABLE_CALIPER OFF CACHE BOOL "" )
set( ENABLE_PAPI OFF CACHE BOOL "" )

set( ENABLE_ESSL OFF CACHE BOOL "" )

set( ENABLE_TRILINOS OFF CACHE BOOL "" )
set( ENABLE_VTK OFF CACHE BOOL "" )

# ROCM options
set( ENABLE_ROCM ON CACHE BOOL "" FORCE )
set( ROCM_ROOT "${HIP_ROOT}" CACHE PATH "" )

# HYPRE options
set( ENABLE_HYPRE_ROCM OFF CACHE BOOL "" )

if( ENABLE_HYPRE_ROCM )
  set( HYPRE_DIR "${GEOSX_TPL_DIR}/hypre-develop-wrtj4ipntfs7a3wvrhp3ibqouije7vcu" CACHE PATH "" )
  set( SUPERLU_DIST_DIR "${GEOSX_TPL_DIR}/superlu-dist-7.2.0-rxdayoia4ct6xim6ki4z2n63w5em3v7h/" CACHE PATH "" ) # TODO: compile superlu@amd with rocm support
else()
  set( HYPRE_DIR "${GEOSX_TPL_DIR}/hypre-develop-26xg3xxx6q3q3ineijpty36eorhsp7mh" CACHE PATH "" )
  set( SUPERLU_DIST_DIR "/sw/crusher/spack-envs/base/opt/cray-sles15-zen3/cce-14.0.0/superlu-dist-7.2.0-kwzr52u5mhumb6exihw2did6h7zkiksq" CACHE PATH "" )
endif()

set( GEOSX_BUILD_OBJ_LIBS OFF CACHE BOOL "" FORCE )

set( ENABLE_GTEST_DEATH_TESTS OFF CACHE BOOL "" )
set( gtest_disable_pthreads ON CACHE BOOL "" )

set( ENABLE_TESTS OFF CACHE BOOL "" FORCE )
set( ENABLE_EXAMPLES OFF CACHE BOOL "" FORCE )
set( ENABLE_BENCHMARKS OFF CACHE BOOL "" FORCE )
set( ENABLE_DOCS OFF CACHE BOOL "" FORCE )