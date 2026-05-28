# Tuolumne toss_4_x86_64_ib_cray llvm-amdgpu@6.4.2 rocm@6.4.2
# MPM Minimal-TPL host-config.
#
# This starts from the generated Tuolumne ROCm/HIP host-config, optionally
# redirects package directories to a flat minimal TPL install, and then applies
# the MPM-only minimal profile.  The Dane CPU minimal host-config is left
# untouched.

include(${CMAKE_CURRENT_LIST_DIR}/tuolumne-toss_4_x86_64_ib_cray-llvm-amdgpu@6.4.2-rocm@6.4.2.cmake)

# Optional minimal-TPL install root.  A flat install/view may be supplied with
# either -DGEOS_TPL_DIR=..., GEOS_MPM_TPL_DIR, GEOS_TPL_DIR, or GEOSX_TPL_DIR.
# Expected subdirectories are:
#   camp raja umpire chai hdf5 conduit silo pugixml fmt vtk
if( NOT DEFINED GEOS_TPL_DIR )
  if( DEFINED ENV{GEOS_MPM_TPL_DIR} )
    set( GEOS_TPL_DIR "$ENV{GEOS_MPM_TPL_DIR}" CACHE PATH "" FORCE )
  elseif( DEFINED ENV{GEOS_TPL_DIR} )
    set( GEOS_TPL_DIR "$ENV{GEOS_TPL_DIR}" CACHE PATH "" FORCE )
  elseif( DEFINED ENV{GEOSX_TPL_DIR} )
    set( GEOS_TPL_DIR "$ENV{GEOSX_TPL_DIR}" CACHE PATH "" FORCE )
  endif()
endif()

function( geos_mpm_set_tpl_dir cache_var subdir )
  if( DEFINED GEOS_TPL_DIR AND EXISTS "${GEOS_TPL_DIR}/${subdir}" )
    set( ${cache_var} "${GEOS_TPL_DIR}/${subdir}" CACHE PATH "" FORCE )
  endif()
endfunction()

geos_mpm_set_tpl_dir( CAMP_DIR camp )
geos_mpm_set_tpl_dir( RAJA_DIR raja )
geos_mpm_set_tpl_dir( UMPIRE_DIR umpire )
geos_mpm_set_tpl_dir( CHAI_DIR chai )
geos_mpm_set_tpl_dir( HDF5_DIR hdf5 )
geos_mpm_set_tpl_dir( CONDUIT_DIR conduit )
geos_mpm_set_tpl_dir( SILO_DIR silo )
geos_mpm_set_tpl_dir( PUGIXML_DIR pugixml )
geos_mpm_set_tpl_dir( FMT_DIR fmt )
geos_mpm_set_tpl_dir( VTK_DIR vtk )

include(${CMAKE_CURRENT_LIST_DIR}/../profiles/mpm-minimal-tpl.cmake)

# Keep the Tuolumne device settings explicit after the minimal profile is applied.
set( ENABLE_HIP ON CACHE BOOL "" FORCE )
set( ENABLE_CUDA OFF CACHE BOOL "" FORCE )
set( CMAKE_HIP_STANDARD "17" CACHE STRING "" FORCE )
# CMake native HIP support rejects the hipcc wrapper as CMAKE_HIP_COMPILER.
# Use ROCm Clang directly and provide the ROCm root for HIP package discovery.
set( CMAKE_HIP_COMPILER "/opt/rocm-6.4.2/llvm/bin/amdclang++" CACHE PATH "" FORCE )
set( CMAKE_HIP_COMPILER_ROCM_ROOT "/opt/rocm-6.4.2" CACHE PATH "" FORCE )
set( CMAKE_HIP_ARCHITECTURES "gfx942" CACHE STRING "" FORCE )
set( ROCM_PATH "/opt/rocm-6.4.2" CACHE PATH "" FORCE )

# Tuolumne jobs are launched with Flux.  Unit/integration tests are disabled by
# the MPM minimal profile, but keep MPIEXEC consistent for any explicit user tests.
set( MPIEXEC "flux" CACHE PATH "" FORCE )
set( MPIEXEC_NUMPROC_FLAG "-n" CACHE STRING "" FORCE )
set( MPIEXEC_PREFLAGS "run;--env=*;--exclusive" CACHE STRING "" FORCE )
