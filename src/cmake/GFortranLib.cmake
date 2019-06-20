##------------------------------------------------------------------------------
##
## Simple script to detect the path for runtime libraries required by gfortran:
## - libgfortran.so
## - libquadmath.so
##
##------------------------------------------------------------------------------
EXECUTE_PROCESS( COMMAND ${CMAKE_C_COMPILER} -print-file-name=libgfortran.so
  OUTPUT_VARIABLE GFORTRANLIB OUTPUT_STRIP_TRAILING_WHITESPACE )
EXECUTE_PROCESS( COMMAND ${CMAKE_C_COMPILER} -print-file-name=libquadmath.so
  OUTPUT_VARIABLE QUADMATHLIB OUTPUT_STRIP_TRAILING_WHITESPACE )

MESSAGE( STATUS "The FORTRAN runtime library are ${GFORTRANLIB} and ${QUADMATHLIB}" )
SET( GEOSX_FORTRAN_RUNTIME ${GFORTRANLIB} ${QUADMATHLIB} CACHE INTERNAL "The fortran runtime library" )
