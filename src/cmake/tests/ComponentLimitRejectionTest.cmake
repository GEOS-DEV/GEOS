# SPDX-License-Identifier: LGPL-2.1-only

cmake_minimum_required( VERSION 3.24 )
execute_process( COMMAND "${GEOS_EXECUTABLE}" -i "${INPUT_FILE}" -o "${OUTPUT_DIRECTORY}" -n componentLimitRejection
                 RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error )
if( result EQUAL 0 OR NOT "${output}${error}" MATCHES "Rebuild with GEOS_MAX_COMPONENTS >= 9" )
  message( FATAL_ERROR "Expected an input error identifying the required component limit.\n${output}${error}" )
endif()
