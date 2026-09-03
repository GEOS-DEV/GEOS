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

set( _geos_sanitizer_uses_asan 0 )
string( JOIN " " _geos_sanitizer_flag_probe
        "${CMAKE_C_FLAGS}" "${CMAKE_CXX_FLAGS}" "${CMAKE_EXE_LINKER_FLAGS}"
        "${CMAKE_SHARED_LINKER_FLAGS}" "${CMAKE_HIP_FLAGS}" "${CMAKE_CUDA_FLAGS}" )
if( _geos_sanitizer_flag_probe MATCHES "-fsanitize=address" )
  set( _geos_sanitizer_uses_asan 1 )
endif()
unset( _geos_sanitizer_flag_probe )

# Defaults match scripts/runIntegratedTests.sh, without log_path so CTest
# captures sanitizer output on the test's stderr.
set( _geos_sanitizer_fast_unwind 0 )
set( _geos_sanitizer_detect_leaks 1 )
if( ENABLE_HIP AND _geos_sanitizer_uses_asan )
  # ROCm's host runtime can allocate during shared-library initialization.
  # ASan's slow unwinder takes libgcc's unwinder lock from that initializer,
  # which can deadlock before GEOS reaches main().
  set( _geos_sanitizer_fast_unwind 1 )
  # LSan cannot complete its stop-the-world scan reliably with the HSA and
  # MPICH helper threads used by ROCm. Avoid hanging every HIP test during
  # the exit-time leak scan when a legacy HIP ASan configuration is supplied.
  set( _geos_sanitizer_detect_leaks 0 )
endif()
set( GEOS_ASAN_OPTIONS
     "protect_shadow_gap=0:abort_on_error=1:detect_leaks=${_geos_sanitizer_detect_leaks}:leak_check_at_exit=${_geos_sanitizer_detect_leaks}:fast_unwind_on_malloc=${_geos_sanitizer_fast_unwind}:print_summary=1:halt_on_error=1"
     CACHE STRING "ASan runtime options used when running GEOS tests" )
set( GEOS_UBSAN_OPTIONS
     "halt_on_error=1:print_stacktrace=1"
     CACHE STRING "UBSan runtime options used when running GEOS tests" )
set( GEOS_LSAN_OPTIONS_BASE
     "detect_leaks=${_geos_sanitizer_detect_leaks}:leak_check_at_exit=${_geos_sanitizer_detect_leaks}:fast_unwind_on_malloc=${_geos_sanitizer_fast_unwind}:print_suppressions=0"
     CACHE STRING "LSan runtime options used when running GEOS tests (suppressions path is appended)" )
set( GEOS_LSAN_OPTIONS
     "${GEOS_LSAN_OPTIONS_BASE}:suppressions=${GEOS_LSAN_SUPPRESSIONS}" )

set( GEOS_SANITIZER_CMAKE_ENV "UBSAN_OPTIONS=${GEOS_UBSAN_OPTIONS}" )
if( _geos_sanitizer_uses_asan )
  list( APPEND GEOS_SANITIZER_CMAKE_ENV
        "ASAN_OPTIONS=${GEOS_ASAN_OPTIONS}"
        "LSAN_OPTIONS=${GEOS_LSAN_OPTIONS}" )
endif()

# Build-time GEOS tools are also used by custom targets such as schema
# generation.  CUDA/HIP host-device compilation can emit duplicate internal
# string globals into one shared library, which ASan reports as an ODR
# violation even though both registrations come from the same object.  Keep
# ODR checking enabled for CTest, but avoid this offload false positive when
# running a build-time tool.
set( GEOS_SANITIZER_BUILD_TOOL_CMAKE_ENV ${GEOS_SANITIZER_CMAKE_ENV} )
if( _geos_sanitizer_uses_asan AND ( ENABLE_CUDA OR ENABLE_HIP ) )
  set( _geos_build_tool_asan_options "${GEOS_ASAN_OPTIONS}:detect_odr_violation=0" )
  set( GEOS_SANITIZER_BUILD_TOOL_CMAKE_ENV
       "ASAN_OPTIONS=${_geos_build_tool_asan_options}"
       "LSAN_OPTIONS=${GEOS_LSAN_OPTIONS}"
       "UBSAN_OPTIONS=${GEOS_UBSAN_OPTIONS}" )
  unset( _geos_build_tool_asan_options )
endif()

message( STATUS "GEOS sanitizer test environment:" )
message( STATUS "  ASan enabled: ${_geos_sanitizer_uses_asan}" )
if( _geos_sanitizer_uses_asan )
  message( STATUS "  ASAN_OPTIONS=${GEOS_ASAN_OPTIONS}" )
  message( STATUS "  LSAN_OPTIONS=${GEOS_LSAN_OPTIONS}" )
endif()
message( STATUS "  UBSAN_OPTIONS=${GEOS_UBSAN_OPTIONS}" )
unset( _geos_sanitizer_uses_asan )
unset( _geos_sanitizer_detect_leaks )
unset( _geos_sanitizer_fast_unwind )

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
