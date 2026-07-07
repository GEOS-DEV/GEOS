# RZHound toss_4_x86_64_ib gcc@12.1.1 MPM Minimal-TPL host-config.
#
# RZHound is a CPU-only TOSS4 Sapphire Rapids Slurm cluster.  Start from the
# Dane CPU host-config wrapper for compiler/MPI options, then let RZHound users
# redirect the GEOS third-party libraries to an RZHound-visible minimal TPL view.
#
# A flat TPL install/view may be supplied with any of:
#   -DGEOS_TPL_DIR=...
#   GEOS_MPM_TPL_DIR=...
#   GEOS_TPL_DIR=...
#   GEOSX_TPL_DIR=...
# Expected subdirectories are:
#   camp raja umpire chai hdf5 conduit silo pugixml fmt vtk
# Direct package overrides via environment variables such as CONDUIT_DIR,
# HDF5_DIR, SILO_DIR, RAJA_DIR, etc. are also honored.

include(${CMAKE_CURRENT_LIST_DIR}/dane-toss_4_x86_64_ib-gcc@12.1.1-mpm-minimal-tpl.cmake)

if( NOT DEFINED GEOS_TPL_DIR OR "${GEOS_TPL_DIR}" STREQUAL "" )
  if( DEFINED ENV{GEOS_MPM_TPL_DIR} )
    set( GEOS_TPL_DIR "$ENV{GEOS_MPM_TPL_DIR}" CACHE PATH "" FORCE )
  elseif( DEFINED ENV{GEOS_TPL_DIR} )
    set( GEOS_TPL_DIR "$ENV{GEOS_TPL_DIR}" CACHE PATH "" FORCE )
  elseif( DEFINED ENV{GEOSX_TPL_DIR} )
    set( GEOS_TPL_DIR "$ENV{GEOSX_TPL_DIR}" CACHE PATH "" FORCE )
  endif()
endif()

function( geos_mpm_rzhound_set_tpl_dir cache_var subdir )
  if( DEFINED GEOS_TPL_DIR AND EXISTS "${GEOS_TPL_DIR}/${subdir}" )
    set( ${cache_var} "${GEOS_TPL_DIR}/${subdir}" CACHE PATH "" FORCE )
  endif()

  if( DEFINED ENV{${cache_var}} AND EXISTS "$ENV{${cache_var}}" )
    set( ${cache_var} "$ENV{${cache_var}}" CACHE PATH "" FORCE )
  endif()
endfunction()

geos_mpm_rzhound_set_tpl_dir( CAMP_DIR camp )
geos_mpm_rzhound_set_tpl_dir( RAJA_DIR raja )
geos_mpm_rzhound_set_tpl_dir( UMPIRE_DIR umpire )
geos_mpm_rzhound_set_tpl_dir( CHAI_DIR chai )
geos_mpm_rzhound_set_tpl_dir( HDF5_DIR hdf5 )
geos_mpm_rzhound_set_tpl_dir( CONDUIT_DIR conduit )
geos_mpm_rzhound_set_tpl_dir( SILO_DIR silo )
geos_mpm_rzhound_set_tpl_dir( PUGIXML_DIR pugixml )
geos_mpm_rzhound_set_tpl_dir( FMT_DIR fmt )
geos_mpm_rzhound_set_tpl_dir( VTK_DIR vtk )
