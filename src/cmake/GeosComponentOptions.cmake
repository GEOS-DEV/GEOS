# SPDX-License-Identifier: LGPL-2.1-only

set( GEOS_MAX_FLUID_COMPONENTS 5 CACHE STRING
     "Maximum number of components instantiated in compositional flow solvers (2-20)" )
if( NOT GEOS_MAX_FLUID_COMPONENTS MATCHES "^([2-9]|1[0-9]|20)$" )
  message( FATAL_ERROR "GEOS_MAX_FLUID_COMPONENTS must be an integer from 2 to 20." )
endif()

# Keep dispatch and explicit instantiations on the same component range.
set( GEOS_COMPONENT_INSTANTIATIONS "" )
foreach( component RANGE 1 ${GEOS_MAX_FLUID_COMPONENTS} )
  string( APPEND GEOS_COMPONENT_INSTANTIATIONS " MACRO( ${component} )" )
endforeach()

# OBL table interpolation uses workspace exponential in the component count.
# Keep its existing upper range while allowing larger EOS fluid models.
set( GEOS_MAX_OBL_COMPONENTS ${GEOS_MAX_FLUID_COMPONENTS} )
if( GEOS_MAX_OBL_COMPONENTS GREATER 9 )
  set( GEOS_MAX_OBL_COMPONENTS 9 )
endif()
set( GEOS_OBL_COMPONENT_INSTANTIATIONS "" )
foreach( component RANGE 1 ${GEOS_MAX_OBL_COMPONENTS} )
  string( APPEND GEOS_OBL_COMPONENT_INSTANTIATIONS " MACRO( ${component} )" )
endforeach()
