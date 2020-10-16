message("\nProcessing GeosxOptions.cmake")
message("CMAKE_HOST_SYSTEM_NAME = ${CMAKE_HOST_SYSTEM_NAME}")
message("CMAKE_SYSTEM_NAME = ${CMAKE_SYSTEM_NAME}")
message("CMAKE_HOST_APPLE = ${CMAKE_HOST_APPLE}")

### OPTIONS ###
option( GEOSX_ENABLE_FPE "" ON)

option( ENABLE_CALIPER "" OFF )

option( ENABLE_MATHPRESSO "" ON )

option( ENABLE_CHAI "Enables CHAI" ON )
option( BUILD_LOCAL_CHAI "Use the local mirrored CHAI" OFF )

option( ENABLE_RAJA "Enables RAJA" ON )
option( BUILD_LOCAL_RAJA "Use the local mirrored RAJA" OFF )
option( RAJA_ENABLE_TBB "" OFF)
option( RAJA_ENABLE_OPENMP "" OFF )
option( RAJA_ENABLE_CUDA "" OFF )
option( RAJA_ENABLE_TESTS "" OFF )

option( ENABLE_GEOSX_PTP "" ON)
option( ENABLE_PAMELA "" ON )
option( ENABLE_PVTPackage "" ON )

option( ENABLE_UNCRUSTIFY "" ON )

option( ENABLE_XML_UPDATES "" ON )

option( ENABLE_FORTRAN "Enables Fortran support" OFF)

option( ENABLE_METIS "Enables METIS" ON )
option( ENABLE_PARMETIS "Enables PARMETIS" ON )

option( ENABLE_VTK "Enables VTK" ON )

option( ENABLE_TOTALVIEW_OUTPUT "Enables Totalview custom view" OFF )

option( ENABLE_SUPERLU_DIST "Enables SUPERLU_DIST" ON )
option( ENABLE_TRILINOS "Enables TRILINOS" ON )
option( ENABLE_HYPRE "Enables HYPRE" ON )
option( ENABLE_PETSC "Enables PETSC" OFF )
option( ENABLE_SUITESPARSE "Enables SUITESPARSE" ON )

#if ( "${CMAKE_HOST_APPLE}" )
#  option( ENABLE_PETSC "Enables PETSC" OFF )
#else()
#  option( ENABLE_PETSC "Enables PETSC" ON )
#endif()

### LAI SETUP ###

set( supported_LAI Trilinos Hypre Petsc )
set( GEOSX_LA_INTERFACE "Trilinos" CACHE STRING "Linear algebra interface to use in solvers" )
message( STATUS "GEOSX_LA_INTERFACE = ${GEOSX_LA_INTERFACE}" )

if( NOT ( GEOSX_LA_INTERFACE IN_LIST supported_LAI ) )
  message( FATAL_ERROR "GEOSX_LA_INTERFACE must be one of: ${supported_LAI}" )
endif()

### MPI/OMP/CUDA SETUP ###

option( ENABLE_MPI "" ON )

option(CUDA_ENABLED "" OFF)

if( CMAKE_HOST_APPLE AND CMAKE_CXX_COMPILER_ID STREQUAL "Clang" )
  option(ENABLE_OPENMP     "Enables OpenMP compiler support" OFF)  
else()
  option(ENABLE_OPENMP     "Enables OpenMP compiler support" ON)
endif()

### BUILD & BLT SETUP ###

option( BUILD_OBJ_LIBS "Builds coreComponent modules as object libraries" OFF)

option( GEOSX_BUILD_SHARED_LIBS "Builds geosx_core as a shared library " ON )

#set(CMAKE_POSITION_INDEPENDENT_CODE ON  CACHE BOOL "" FORCE)
#blt_append_custom_compiler_flag(FLAGS_VAR CMAKE_CXX_FLAGS DEFAULT -rdynamic)
#set(CMAKE_EXE_LINKER_FLAGS "-rdynamic")

#set(SPHINX_EXECUTABLE "sphinx-build" CACHE PATH "")
#include(cmake/blt/cmake/thirdparty/FindSphinx.cmake)
#message( "SPHINX_FOUND = ${SPHINX_FOUND}" )
#message( "SPHINX_EXECUTABLE = ${SPHINX_EXECUTABLE}" )

if(NOT BLT_CXX_STD STREQUAL c++14)
    MESSAGE(FATAL_ERROR "c++14 is NOT enabled. GEOSX requires c++14")
endif(NOT BLT_CXX_STD STREQUAL c++14)

message("CMAKE_CXX_COMPILER_ID = ${CMAKE_CXX_COMPILER_ID}")

blt_append_custom_compiler_flag( FLAGS_VAR CMAKE_CXX_FLAGS DEFAULT "${OpenMP_CXX_FLAGS}")
blt_append_custom_compiler_flag( FLAGS_VAR CMAKE_CXX_FLAGS
                                 GNU   "-Wall -Wextra -Wpedantic -pedantic-errors -Wshadow -Wfloat-equal -Wcast-align -Wcast-qual"
                                 CLANG "-Wall -Wextra -Wpedantic -pedantic-errors -Wshadow -Wfloat-equal -Wcast-align -Wcast-qual"
                               )

blt_append_custom_compiler_flag( FLAGS_VAR CMAKE_CXX_FLAGS_DEBUG
                                 GNU ""
                                 CLANG "-fstandalone-debug"
                                )

blt_append_custom_compiler_flag(FLAGS_VAR GEOSX_NINJA_FLAGS
                  DEFAULT     " "
                  GNU         "-fdiagnostics-color=always"
                  CLANG       "-fcolor-diagnostics"
                  )

if( ${CMAKE_MAKE_PROGRAM} STREQUAL "ninja" OR ${CMAKE_MAKE_PROGRAM} MATCHES ".*/ninja$" )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GEOSX_NINJA_FLAGS}")
endif()

if( CMAKE_HOST_APPLE )
    set(GEOSX_LINK_PREPEND_FLAG "-Wl,-force_load" CACHE STRING "")
    set(GEOSX_LINK_POSTPEND_FLAG "" CACHE STRING "")
elseif( CUDA_ENABLED )
    set(GEOSX_LINK_PREPEND_FLAG  "-Xcompiler \\\\\"-Wl,--whole-archive\\\\\""    CACHE STRING "")
    set(GEOSX_LINK_POSTPEND_FLAG "-Xcompiler \\\\\"-Wl,--no-whole-archive\\\\\"" CACHE STRING "")
else()
    set(GEOSX_LINK_PREPEND_FLAG  "-Wl,--whole-archive"    CACHE STRING "")
    set(GEOSX_LINK_POSTPEND_FLAG "-Wl,--no-whole-archive" CACHE STRING "")
endif()

message("CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")
message( "GEOSX_LINK_PREPEND_FLAG=${GEOSX_LINK_PREPEND_FLAG}" )
message( "GEOSX_LINK_POSTPEND_FLAG=${GEOSX_LINK_POSTPEND_FLAG}" )
message("Leaving GeosxOptions.cmake\n")
