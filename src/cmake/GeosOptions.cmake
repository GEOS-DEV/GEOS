message( "\nProcessing GeosOptions.cmake" )
message( "CMAKE_HOST_SYSTEM_NAME = ${CMAKE_HOST_SYSTEM_NAME}" )
message( "CMAKE_SYSTEM_NAME = ${CMAKE_SYSTEM_NAME}" )
message( "CMAKE_HOST_APPLE = ${CMAKE_HOST_APPLE}" )

### OPTIONS ###
option( GEOS_ENABLE_FPE "Enables floating point exceptions" ON )
option( GEOS_ENABLE_TESTS "Enables unit tests" ON )
option( GEOS_ENABLE_SANITIZERS "Build and test with AddressSanitizer and UndefinedBehaviorSanitizer" OFF )
option( ENABLE_CALIPER "Enables Caliper instrumentation" OFF )

option( ENABLE_MATHPRESSO "" ON )

option( ENABLE_CHAI "Enables CHAI" ON )
option( BUILD_LOCAL_CHAI "Use the local mirrored CHAI" OFF )

option( ENABLE_RAJA "Enables RAJA" ON )
option( BUILD_LOCAL_RAJA "Use the local mirrored RAJA" OFF )
option( RAJA_ENABLE_TBB "" OFF )
option( RAJA_ENABLE_OPENMP "" OFF )
option( RAJA_ENABLE_CUDA "" OFF )
option( RAJA_ENABLE_HIP "" OFF )
option( RAJA_ENABLE_TESTS "" OFF )

option( GEOS_USE_GIT_VERSION_INFO "Enables querying git for version metadata" ON )

option( GEOS_ENABLE_BOUNDS_CHECK "Enables array bounds checking" OFF )
if( NOT CMAKE_CONFIGURATION_TYPES )
    ######################################################
    # Add define we can use when debug builds are enabled
    ######################################################
    if ( CMAKE_BUILD_TYPE MATCHES "Debug" )
        set( GEOS_ENABLE_BOUNDS_CHECK ON CACHE BOOL "" FORCE )
    endif()
endif()
if( GEOS_ENABLE_BOUNDS_CHECK )
  set( LVARRAY_BOUNDS_CHECK ON CACHE BOOL "" FORCE )
endif()

option( ENABLE_HPCREACT "" ON )

option( ENABLE_UNCRUSTIFY "" ON )

option( ENABLE_XML_UPDATES "" ON )

option( ENABLE_FORTRAN "Enables Fortran support" OFF )

option( ENABLE_METIS "Enables METIS" ON )
option( ENABLE_PARMETIS "Enables PARMETIS" ON )
option( ENABLE_SCOTCH "Enables SCOTCH" ON )

option( ENABLE_SILO "Enables SILO output" ON )
option( ENABLE_VTK "Enables VTK" ON )

option( ENABLE_TOTALVIEW_OUTPUT "Enables Totalview custom view" OFF )

option( ENABLE_SUPERLU_DIST "Enables SUPERLU_DIST" ON )
option( ENABLE_TRILINOS "Enables TRILINOS" ON )
option( ENABLE_HYPRE "Enables HYPRE" ON )
option( ENABLE_HYPREDRV "Enables HYPREDRV" OFF )
option( ENABLE_PETSC "Enables PETSC" OFF )
option( ENABLE_SUITESPARSE "Enables SUITESPARSE" ON )

set( HYPREDRV_DIR "" CACHE PATH "Path to a HYPREDRV installation prefix or package config directory" )

option( ENABLE_HYPRE_MIXINT "Enables mixed int32/int64 local/global" ON )

set( HYPRE_DEVICE_OPTIONS CPU CUDA HIP )
if( NOT ENABLE_HYPRE_DEVICE )
  set( ENABLE_HYPRE_DEVICE CPU )
endif()
if(NOT ${ENABLE_HYPRE_DEVICE} IN_LIST HYPRE_DEVICE_OPTIONS )
    message(FATAL_ERROR "Set ENABLE_HYPRE_DEVICE to CPU, CUDA, or HIP.")
endif()

#if ( "${CMAKE_HOST_APPLE}" )
#  option( ENABLE_PETSC "Enables PETSC" OFF )
#else()
#  option( ENABLE_PETSC "Enables PETSC" ON )
#endif()

### LAI SETUP ###

set( supported_LAI Trilinos Hypre Petsc )
set( GEOS_LA_INTERFACE "Hypre" CACHE STRING "Linear algebra interface to use in solvers" )
message( STATUS "GEOS_LA_INTERFACE = ${GEOS_LA_INTERFACE}" )

if( NOT ( GEOS_LA_INTERFACE IN_LIST supported_LAI ) )
  message( FATAL_ERROR "GEOS_LA_INTERFACE must be one of: ${supported_LAI}" )
endif()

### MPI/OMP/CUDA/HIP SETUP ###

option( ENABLE_MPI "" ON )

option( ENABLE_CUDA "" OFF )

option( ENABLE_HIP "" OFF )

if( CMAKE_HOST_APPLE AND CMAKE_CXX_COMPILER_ID STREQUAL "Clang" )
  option( ENABLE_OPENMP "Enables OpenMP compiler support" OFF )
else()
  option( ENABLE_OPENMP "Enables OpenMP compiler support" ON )
endif()

option( ENABLE_CUDA_STACK_SIZE "Allows the CUDA stack size limit to be adjusted" OFF )

### BUILD & BLT SETUP ###

option( GEOS_INSTALL_SCHEMA "Enables schema generation and installation" ON )

option( GEOS_BUILD_OBJ_LIBS "Builds coreComponent modules as object libraries" OFF )

option( GEOS_BUILD_SHARED_LIBS "Builds geosx_core as a shared library " ON )

set( GEOS_PARALLEL_COMPILE_JOBS "" CACHE STRING "Maximum number of concurrent compilation jobs" )
if( GEOS_PARALLEL_COMPILE_JOBS )
    set_property( GLOBAL APPEND PROPERTY JOB_POOLS compile_job_pool=${GEOS_PARALLEL_COMPILE_JOBS} )
    set( CMAKE_JOB_POOL_COMPILE compile_job_pool )
endif()

set( GEOS_PARALLEL_LINK_JOBS "" CACHE STRING "Maximum number of concurrent link jobs" )
if( GEOS_PARALLEL_LINK_JOBS )
    set_property( GLOBAL APPEND PROPERTY JOB_POOLS link_job_pool=${GEOS_PARALLEL_LINK_JOBS} )
    set( CMAKE_JOB_POOL_LINK link_job_pool )
endif()

# Physics packages
option( GEOS_ENABLE_CONTACT "Enables contact physics package" ON )
option( GEOS_ENABLE_FLUIDFLOW "Enables fluid flow physics package" ON )
option( GEOS_ENABLE_INDUCEDSEISMICITY "Enables induced seismicity physics package" ON )
option( GEOS_ENABLE_MULTIPHYSICS "Enables multiphysics physics package" ON )
option( GEOS_ENABLE_SIMPLEPDE "Enables simple PDE physics package" ON )
option( GEOS_ENABLE_SOLIDMECHANICS "Enables solid mechanics physics package" ON )
option( GEOS_ENABLE_SURFACEGENERATION "Enables surface generation physics package" ON )
option( GEOS_ENABLE_WAVEPROPAGATION "Enables wave propagation physics package" ON )

#set(CMAKE_POSITION_INDEPENDENT_CODE ON  CACHE BOOL "" FORCE)
#blt_append_custom_compiler_flag(FLAGS_VAR CMAKE_CXX_FLAGS DEFAULT -rdynamic)
#set(CMAKE_EXE_LINKER_FLAGS "-rdynamic")

#set(SPHINX_EXECUTABLE "sphinx-build" CACHE PATH "")
#include(cmake/blt/cmake/thirdparty/FindSphinx.cmake)
#message( "SPHINX_FOUND = ${SPHINX_FOUND}" )
#message( "SPHINX_EXECUTABLE = ${SPHINX_EXECUTABLE}" )

if( NOT BLT_CXX_STD STREQUAL c++20 )
    MESSAGE( FATAL_ERROR "c++20 is NOT enabled. GEOS requires c++20" )
endif( NOT BLT_CXX_STD STREQUAL c++20 )

message( "CMAKE_CXX_COMPILER_ID = ${CMAKE_CXX_COMPILER_ID}" )

blt_append_custom_compiler_flag( FLAGS_VAR CMAKE_CXX_FLAGS DEFAULT "${OpenMP_CXX_FLAGS}" )
blt_append_custom_compiler_flag( FLAGS_VAR CMAKE_CXX_FLAGS
                                 GNU   "-Wpedantic -pedantic-errors -Wshadow -Wfloat-equal -Wcast-align -Wcast-qual"
                                 CLANG "-Wpedantic -pedantic-errors -Wshadow -Wfloat-equal -Wno-cast-align -Wcast-qual"
                               )

if( ENABLE_HIP )
  # amdclang compiles C++ sources through the HIP driver.  GEOS has existing
  # [=] lambdas that implicitly capture this; C++20 diagnoses that pattern.
  blt_append_custom_compiler_flag( FLAGS_VAR CMAKE_CXX_FLAGS
                                   CLANG "-Wno-deprecated-this-capture -Wno-unused-parameter -Wno-unused-variable -Wno-unused-lambda-capture -Wno-gpu-maybe-wrong-side" )
endif()

blt_append_custom_compiler_flag( FLAGS_VAR CMAKE_CXX_FLAGS_DEBUG
                                 GNU "-Wno-unused-parameter -Wno-unused-variable"
                                 CLANG "-Wno-unused-parameter -Wno-unused-variable -fstandalone-debug"
                               )

blt_append_custom_compiler_flag( FLAGS_VAR GEOS_NINJA_FLAGS
                                 DEFAULT " "
                                 GNU     "-fdiagnostics-color=always"
                                 CLANG   "-fcolor-diagnostics"
                               )

# clang-13 and gcc complains about unused-but-set variable.
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-Wunused-but-set-variable" CXX_UNUSED_BUT_SET_VAR)
# clang-22 reports google-benchmark's use of __COUNTER__ inside a #if directive as a
# C2y extension. The check keeps the flag off compilers that do not know it, where
# an unrecognized -Wno-* would itself become an error under ENABLE_WARNINGS_AS_ERRORS.
CHECK_CXX_COMPILER_FLAG("-Wc2y-extensions" CXX_C2Y_EXTENSIONS)
if (ENABLE_GBENCHMARK)
    blt_add_target_compile_flags(TO benchmark
                                FLAGS $<$<AND:$<BOOL:${CXX_UNUSED_BUT_SET_VAR}>,$<COMPILE_LANGUAGE:CXX>>:-Wno-unused-but-set-variable>
                                      $<$<AND:$<BOOL:${CXX_C2Y_EXTENSIONS}>,$<COMPILE_LANGUAGE:CXX>>:-Wno-c2y-extensions>
                                )
endif()

if( GEOS_ENABLE_FPE )
  # amdclang compiles HIP sources through the CXX driver.  The host-only
  # floating-point exception flag is rejected for the device compilation.
  if( ENABLE_HIP )
    message( STATUS "GEOS_ENABLE_FPE is ON, but HIP builds do not support -ffp-exception-behavior=strict; skipping the flag." )
  else()
    check_cxx_compiler_flag( "-ffp-exception-behavior=strict" GEOS_CXX_HAS_FP_EXCEPTION_BEHAVIOR_STRICT)
    if( GEOS_CXX_HAS_FP_EXCEPTION_BEHAVIOR_STRICT )
      blt_append_custom_compiler_flag( FLAGS_VAR CMAKE_CXX_FLAGS CLANG "-ffp-exception-behavior=strict" )
    else()
      message( WARNING "GEOS_ENABLE_FPE is ON, but ${CMAKE_CXX_COMPILER_ID} does not support -ffp-exception-behavior=strict." )
    endif()
  endif()
endif()

if( ${CMAKE_MAKE_PROGRAM} STREQUAL "ninja" OR ${CMAKE_MAKE_PROGRAM} MATCHES ".*/ninja$" )
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GEOS_NINJA_FLAGS}" )
endif()

set( GEOS_LOCALINDEX_TYPE "int" CACHE STRING "" )
if( ENABLE_HYPRE_MIXINT )
  set( GEOS_GLOBALINDEX_TYPE "long long int" CACHE STRING "" )
else()
  set( GEOS_GLOBALINDEX_TYPE "int" CACHE STRING "" )
endif()

if( GEOS_LOCALINDEX_TYPE STREQUAL "int" )
    set( GEOS_LOCALINDEX_TYPE_FLAG "0" CACHE STRING "" FORCE )
elseif( GEOS_LOCALINDEX_TYPE STREQUAL "long int" )
    set( GEOS_LOCALINDEX_TYPE_FLAG "1" CACHE STRING "" FORCE )
elseif( GEOS_LOCALINDEX_TYPE STREQUAL "long long int" )
    set( GEOS_LOCALINDEX_TYPE_FLAG "2" CACHE STRING "" FORCE )
elseif( GEOS_LOCALINDEX_TYPE STREQUAL "std::ptrdiff_t" )
    set( GEOS_LOCALINDEX_TYPE_FLAG "3" CACHE STRING "" FORCE )
else( TRUE )
    message( FATAL_ERROR "GEOS_LOCALINDEX_TYPE_FLAG not set for ${GEOS_LOCALINDEX_TYPE}" )
endif()



if( GEOS_GLOBALINDEX_TYPE STREQUAL "int" )
    set( GEOS_GLOBALINDEX_TYPE_FLAG "0" CACHE STRING "" FORCE )
elseif( GEOS_GLOBALINDEX_TYPE STREQUAL "long int" )
    set( GEOS_GLOBALINDEX_TYPE_FLAG "1" CACHE STRING "" FORCE )
elseif( GEOS_GLOBALINDEX_TYPE STREQUAL "long long int" )
    set( GEOS_GLOBALINDEX_TYPE_FLAG "2" CACHE STRING "" FORCE )
else( TRUE )
    message( FATAL_ERROR "GEOS_GLOBALINDEX_TYPE_FLAG not set for ${GEOS_GLOBALINDEX_TYPE}" )
endif()

set( GEOS_BLOCK_SIZE 32 )
if( ENABLE_CUDA )
  set( GEOS_BLOCK_SIZE 32 )
endif()
if( ENABLE_HIP )
  set( GEOS_BLOCK_SIZE 64 )
endif()

message( "localIndex is an alias for ${GEOS_LOCALINDEX_TYPE}" )
message( "globalIndex is an alias for ${GEOS_GLOBALINDEX_TYPE}" )
message( "GEOS_LOCALINDEX_TYPE_FLAG = ${GEOS_LOCALINDEX_TYPE_FLAG}" )
message( "GEOS_GLOBALINDEX_TYPE_FLAG = ${GEOS_GLOBALINDEX_TYPE_FLAG}" )


message( "CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}" )

# A TPL host-config or -DCMAKE_*_FLAGS=-fsanitize=... already instruments the
# build. Treat that as a sanitizer build even if the cache option was left off.
if( NOT GEOS_ENABLE_SANITIZERS )
  string( JOIN " " _geos_sanitizer_flag_probe
          "${CMAKE_C_FLAGS}" "${CMAKE_CXX_FLAGS}"
          "${CMAKE_EXE_LINKER_FLAGS}" "${CMAKE_SHARED_LINKER_FLAGS}" )
  if( _geos_sanitizer_flag_probe MATCHES "-fsanitize=" )
    set( GEOS_ENABLE_SANITIZERS ON CACHE BOOL
         "Detected from compiler/linker flags" FORCE )
    message( STATUS "GEOS_ENABLE_SANITIZERS enabled because sanitizer flags are present" )
  endif()
  unset( _geos_sanitizer_flag_probe )
endif()

if( GEOS_ENABLE_SANITIZERS AND ENABLE_HIP AND CMAKE_CXX_COMPILER_ID STREQUAL "Clang"
    AND CMAKE_CXX_FLAGS MATCHES "-fsanitize=address" )
  # amdclang emits duplicate ASan metadata for local string constants when
  # host code is compiled as part of a HIP build.  The resulting ODR report
  # is a compiler/runtime false positive and aborts every test before main().
  # Keep ASan's stack/heap instrumentation and UBSan checks while disabling
  # only global instrumentation for the HIP host compiler.
  check_cxx_compiler_flag( "-mllvm=-asan-globals=0" GEOS_CXX_HAS_NO_ASAN_GLOBALS )
  if( GEOS_CXX_HAS_NO_ASAN_GLOBALS )
    # Keep the LLVM backend option off link commands.  This also migrates
    # build directories created with the earlier flag-based workaround.
    string( REPLACE "-mllvm=-asan-globals=0" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}" )
    string( REPLACE "-fno-sanitize-address-use-odr-indicator" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}" )
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}" CACHE STRING "" FORCE )
    add_compile_options( $<$<COMPILE_LANG_AND_ID:CXX,Clang>:-mllvm=-asan-globals=0> )

    if( CMAKE_HIP_FLAGS MATCHES "-fsanitize=address" AND NOT CMAKE_HIP_FLAGS MATCHES "-asan-globals=0" )
      set( CMAKE_HIP_FLAGS "${CMAKE_HIP_FLAGS} -Xarch_host -mllvm=-asan-globals=0" )
      set( CMAKE_HIP_FLAGS "${CMAKE_HIP_FLAGS}" CACHE STRING "" FORCE )
    endif()
    unset( _geos_cxx_no_asan_globals_pos )
  endif()
endif()

if( GEOS_ENABLE_SANITIZERS AND ENABLE_VTK )
  # VTK and GEOS' VTK mesh registrations construct polymorphic objects during
  # shared-library startup. UBSan's vptr check diagnoses that valid cross-DSO
  # initialization before main(), so disable only vptr instrumentation for
  # the VTK-enabled sanitizer build while retaining ASan and the other UBSan
  # checks.
  set( _geos_vtk_ubsan_flag "-fno-sanitize=vptr" )
  string( FIND "${CMAKE_CXX_FLAGS}" "${_geos_vtk_ubsan_flag}" _geos_vtk_ubsan_pos )
  if( _geos_vtk_ubsan_pos EQUAL -1 )
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${_geos_vtk_ubsan_flag}" )
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}" CACHE STRING "" FORCE )
  endif()
  unset( _geos_vtk_ubsan_flag )
  unset( _geos_vtk_ubsan_pos )
endif()

message( "Leaving GeosOptions.cmake\n" )
