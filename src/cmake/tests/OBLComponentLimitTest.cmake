# SPDX-License-Identifier: LGPL-2.1-only

cmake_minimum_required( VERSION 3.24 )
file( READ "${INPUT_FILE}" input )
get_filename_component( inputDirectory "${INPUT_FILE}" DIRECTORY )
string( REPLACE "./deadoil_3ph_staircase_obl_3d_base.xml"
                "${inputDirectory}/deadoil_3ph_staircase_obl_3d_base.xml" input "${input}" )
string( REPLACE "obl_do_static.txt" "${inputDirectory}/obl_do_static.txt" input "${input}" )
string( REPLACE "numComponents=\"3\"" "numComponents=\"10\"" input "${input}" )
set( rejectionInput "${OUTPUT_DIRECTORY}/oblComponentLimit.xml" )
file( WRITE "${rejectionInput}" "${input}" )
execute_process( COMMAND "${GEOS_EXECUTABLE}" -i "${rejectionInput}" -o "${OUTPUT_DIRECTORY}" -n oblComponentLimit
                 RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error )
if( result EQUAL 0 OR NOT "${output}${error}" MATCHES "OBL table interpolation supports at most nine components" )
  message( FATAL_ERROR "Expected an input error identifying the OBL component limit.\n${output}${error}" )
endif()
