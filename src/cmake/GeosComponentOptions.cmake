# SPDX-License-Identifier: LGPL-2.1-only

set( GEOS_MAX_COMPONENTS 5 CACHE STRING
     "Maximum number of components instantiated in compositional flow solvers (2-9)" )
if( NOT GEOS_MAX_COMPONENTS MATCHES "^[2-9]$" )
  message( FATAL_ERROR "GEOS_MAX_COMPONENTS must be an integer from 2 to 9." )
endif()

# Keep dispatch and explicit instantiations on the same component range.
set( GEOS_COMPONENT_INSTANTIATIONS "" )
foreach( component RANGE 1 ${GEOS_MAX_COMPONENTS} )
  string( APPEND GEOS_COMPONENT_INSTANTIATIONS " MACRO( ${component} )" )
endforeach()
