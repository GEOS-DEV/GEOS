# SPDX-License-Identifier: LGPL-2.1-only

set( GEOS_MAX_FLUID_COMPONENTS 5 CACHE STRING
     "Maximum number of components instantiated in compositional flow solvers (3-20)" )
# The component visitor includes black-oil fluids, which require three
# components. Smaller fluid models remain available with this build limit.
if( NOT GEOS_MAX_FLUID_COMPONENTS MATCHES "^([3-9]|1[0-9]|20)$" )
  message( FATAL_ERROR "GEOS_MAX_FLUID_COMPONENTS must be an integer from 3 to 20." )
endif()

# Keep dispatch and explicit instantiations on the same component range.
set( GEOS_COMPONENT_INSTANTIATIONS "" )
foreach( component RANGE 1 ${GEOS_MAX_FLUID_COMPONENTS} )
  string( APPEND GEOS_COMPONENT_INSTANTIATIONS " MACRO( ${component} )" )
endforeach()

# OBL table interpolation keeps a per-thread workspace of (2^(NC+1)-1) x numOps
# doubles. Seven components with energy need about 307 KiB, below the 512 KiB
# CUDA local-memory limit; eight would need about 679 KiB.
set( GEOS_OBL_COMPONENT_CAP 7 )
set( GEOS_MAX_OBL_COMPONENTS ${GEOS_MAX_FLUID_COMPONENTS} )
if( GEOS_MAX_OBL_COMPONENTS GREATER GEOS_OBL_COMPONENT_CAP )
  set( GEOS_MAX_OBL_COMPONENTS ${GEOS_OBL_COMPONENT_CAP} )
endif()
set( GEOS_OBL_COMPONENT_INSTANTIATIONS "" )
foreach( component RANGE 1 ${GEOS_MAX_OBL_COMPONENTS} )
  string( APPEND GEOS_OBL_COMPONENT_INSTANTIATIONS " MACRO( ${component} )" )
endforeach()
