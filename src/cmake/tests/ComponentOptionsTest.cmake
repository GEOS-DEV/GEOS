# SPDX-License-Identifier: LGPL-2.1-only

cmake_minimum_required( VERSION 3.24 )

set( optionsFile "${CMAKE_CURRENT_LIST_DIR}/../GeosComponentOptions.cmake" )
include( "${optionsFile}" )
if( NOT GEOS_MAX_COMPONENTS EQUAL 5 )
  message( FATAL_ERROR "The default component limit must remain five." )
endif()

foreach( limit RANGE 2 9 )
  set( GEOS_MAX_COMPONENTS ${limit} CACHE STRING "" FORCE )
  include( "${optionsFile}" )
  string( REGEX MATCHALL "MACRO\\( [0-9]+ \\)" counts "${GEOS_COMPONENT_INSTANTIATIONS}" )
  list( LENGTH counts count )
  if( NOT count EQUAL limit )
    message( FATAL_ERROR "Wrong number of instantiated component counts for limit ${limit}." )
  endif()
  foreach( component RANGE 1 ${limit} )
    if( NOT "MACRO( ${component} )" IN_LIST counts )
      message( FATAL_ERROR "Missing component count ${component} for limit ${limit}." )
    endif()
  endforeach()
endforeach()

foreach( invalid IN ITEMS 0 1 10 -1 2.5 abc 02 "" )
  execute_process( COMMAND "${CMAKE_COMMAND}" "-DGEOS_MAX_COMPONENTS=${invalid}" -P "${optionsFile}"
                   RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error )
  if( result EQUAL 0 OR NOT error MATCHES "GEOS_MAX_COMPONENTS must be an integer from 2 to 9" )
    message( FATAL_ERROR "Invalid component limit '${invalid}' was not rejected: ${output}${error}" )
  endif()
endforeach()
