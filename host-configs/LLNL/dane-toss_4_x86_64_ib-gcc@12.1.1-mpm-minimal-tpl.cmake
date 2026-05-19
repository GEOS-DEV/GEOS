# Dane toss_4_x86_64_ib gcc@12.1.1 MPM Minimal-TPL host-config.
# This starts from the generated Dane host-config and then strips non-MPM packages/TPLs.

include(${CMAKE_CURRENT_LIST_DIR}/dane-toss_4_x86_64_ib-gcc@12.1.1.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/../profiles/mpm-minimal-tpl.cmake)

# This is the CPU-only Dane MPM configuration. Force device backends off so
# a stale cache or Tuolumne environment cannot make the Dane build create BLT
# HIP/CUDA targets.
set( ENABLE_HIP OFF CACHE BOOL "" FORCE )
set( ENABLE_CUDA OFF CACHE BOOL "" FORCE )
unset( CMAKE_HIP_COMPILER CACHE )
unset( CMAKE_HIP_COMPILER_ROCM_ROOT CACHE )
unset( CMAKE_HIP_ARCHITECTURES CACHE )
unset( ROCM_PATH CACHE )
