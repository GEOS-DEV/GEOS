# Runtime environment for ASan/LSan/UBSan when running GEOS tests.
# Applied to every CTest test and to the ATS helper targets.

if( NOT GEOS_ENABLE_SANITIZERS )
  function( geos_apply_sanitizer_test_environment )
  endfunction()
  set( GEOS_SANITIZER_CMAKE_ENV )
  return()
endif()

set( GEOS_LSAN_SUPPRESSIONS "${CMAKE_BINARY_DIR}/lsan.supp" )
set( _geos_lsan_supp_src "${CMAKE_SOURCE_DIR}/../scripts/lsan.supp" )
if( NOT EXISTS "${_geos_lsan_supp_src}" )
  message( FATAL_ERROR "GEOS_ENABLE_SANITIZERS is ON but ${_geos_lsan_supp_src} was not found" )
endif()
configure_file( "${_geos_lsan_supp_src}" "${GEOS_LSAN_SUPPRESSIONS}" COPYONLY )
unset( _geos_lsan_supp_src )

# Defaults match scripts/runIntegratedTests.sh, without log_path so CTest
# captures sanitizer output on the test's stderr.
set( GEOS_ASAN_OPTIONS
     "abort_on_error=1:detect_leaks=1:leak_check_at_exit=1:fast_unwind_on_malloc=0:print_summary=1:halt_on_error=1"
     CACHE STRING "ASan runtime options used when running GEOS tests" )
set( GEOS_UBSAN_OPTIONS
     "halt_on_error=1:print_stacktrace=1"
     CACHE STRING "UBSan runtime options used when running GEOS tests" )
set( GEOS_LSAN_OPTIONS_BASE
     "detect_leaks=1:leak_check_at_exit=1:fast_unwind_on_malloc=0:print_suppressions=0"
     CACHE STRING "LSan runtime options used when running GEOS tests (suppressions path is appended)" )
set( GEOS_LSAN_OPTIONS
     "${GEOS_LSAN_OPTIONS_BASE}:suppressions=${GEOS_LSAN_SUPPRESSIONS}" )

set( GEOS_SANITIZER_CMAKE_ENV
     "ASAN_OPTIONS=${GEOS_ASAN_OPTIONS}"
     "LSAN_OPTIONS=${GEOS_LSAN_OPTIONS}"
     "UBSAN_OPTIONS=${GEOS_UBSAN_OPTIONS}" )

message( STATUS "GEOS sanitizer test environment:" )
message( STATUS "  ASAN_OPTIONS=${GEOS_ASAN_OPTIONS}" )
message( STATUS "  LSAN_OPTIONS=${GEOS_LSAN_OPTIONS}" )
message( STATUS "  UBSAN_OPTIONS=${GEOS_UBSAN_OPTIONS}" )

function( geos_apply_sanitizer_test_environment )
  if( NOT GEOS_ENABLE_SANITIZERS )
    return()
  endif()
  geos_apply_sanitizer_test_environment_dir( "${CMAKE_SOURCE_DIR}" )
endfunction()

function( geos_apply_sanitizer_test_environment_dir dir )
  get_property( _geos_sanitizer_tests DIRECTORY "${dir}" PROPERTY TESTS )
  if( _geos_sanitizer_tests )
    set_property( TEST ${_geos_sanitizer_tests} DIRECTORY "${dir}"
                  APPEND PROPERTY ENVIRONMENT ${GEOS_SANITIZER_CMAKE_ENV} )
  endif()
  get_property( _geos_sanitizer_subdirs DIRECTORY "${dir}" PROPERTY SUBDIRECTORIES )
  foreach( _geos_sanitizer_subdir IN LISTS _geos_sanitizer_subdirs )
    geos_apply_sanitizer_test_environment_dir( "${_geos_sanitizer_subdir}" )
  endforeach()
endfunction()
