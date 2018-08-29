message("\nProcessing geosxOptions.cmake")

message("CMAKE_HOST_SYSTEM_NAME = ${CMAKE_HOST_SYSTEM_NAME}")
message("CMAKE_SYSTEM_NAME = ${CMAKE_SYSTEM_NAME}")
message("CMAKE_HOST_APPLE = ${CMAKE_HOST_APPLE}")





# OPTIONS

option( ENABLE_CALIPER "" OFF )
set( CALIPER_INSTALL "" CACHE PATH "")

option( ENABLE_MATHPRESSO "" ON )

option( ENABLE_CHAI "Enables CHAI" ON )
option( BUILD_LOCAL_CHAI "Use the local mirrored CHAI" OFF )

option( ENABLE_RAJA "Enables RAJA" ON )
option( BUILD_LOCAL_RAJA "Use the local mirrored RAJA" OFF )
option( RAJA_ENABLE_TBB "" OFF)
option( RAJA_ENABLE_OPENMP "" OFF )
option( RAJA_ENABLE_CUDA "" OFF )
option( RAJA_ENABLE_TESTS "" OFF )


option( ENABLE_FPARSER "Enables FPARSER" OFF )

option( ENABLE_UNCRUSTIFY "" ON )

option( ENABLE_FORTRAN "Enables Fortran support" OFF)


option(ENABLE_CONTAINERARRAY_RETURN_PTR     "Enables ViewWrapper to return pointers instead of references" ON )

option( ENABLE_HYPRE "Enables HYPRE" OFF )

option( ENABLE_MPI "" ON )

option(CUDA_ENABLED "" OFF)

if( CMAKE_HOST_APPLE AND CMAKE_CXX_COMPILER_ID STREQUAL "Clang" )
  option(ENABLE_OPENMP     "Enables OpenMP compiler support" OFF)  
else()
  option(ENABLE_OPENMP     "Enables OpenMP compiler support" ON)
endif()


#set( BUILD_SHARED_LIBS ON CACHE PATH "" FORCE)
#set( ENABLE_SHARED_LIBS ON CACHE PATH "" FORCE )

#set(CMAKE_POSITION_INDEPENDENT_CODE ON  CACHE BOOL "" FORCE)
#blt_append_custom_compiler_flag(FLAGS_VAR CMAKE_CXX_FLAGS DEFAULT -rdynamic)
#set(CMAKE_EXE_LINKER_FLAGS "-rdynamic")



set( GEOSX_TPL_ROOT_DIR "../../thirdPartyLibs/" CACHE PATH "" )

message( "GEOSX_TPL_ROOT_DIR = ${GEOSX_TPL_ROOT_DIR}" )
get_filename_component(GEOSX_TPL_ROOT_DIR ${GEOSX_TPL_ROOT_DIR} ABSOLUTE)
message( "Absolute GEOSX_TPL_ROOT_DIR = ${GEOSX_TPL_ROOT_DIR}" )


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
                                 GNU "-Wall -pedantic-errors -Wno-abi -Wextra  -Wshadow -Wfloat-equal	-Wcast-align	-Wpointer-arith	-Wwrite-strings	-Wcast-qual	-Wswitch-default  -Wno-vla  -Wno-switch-default  -Wno-unused-parameter  -Wno-unused-variable  -Wno-unused-function" 
                                 CLANG "-Weverything -Wno-c++98-compat -Wno-c++98-compat-pedantic -Wno-padded -Wno-missing-prototypes -Wno-covered-switch-default -Wno-double-promotion -Wno-documentation -Wno-switch-enum -Wno-sign-conversion -Wno-unused-parameter -Wno-unused-variable -Wno-reserved-id-macro -Wno-weak-vtables -Wno-undefined-func-template" 
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
    set(GEOSX_LINK_PREPEND_FLAG "-Wl,-force_load" CACHE PATH "" FORCE)
    set(GEOSX_LINK_POSTPEND_FLAG "" CACHE PATH "" FORCE)
else()
    set(GEOSX_LINK_PREPEND_FLAG  "-Wl,--whole-archive"    CACHE PATH "" FORCE)
    set(GEOSX_LINK_POSTPEND_FLAG "-Wl,--no-whole-archive" CACHE PATH "" FORCE)
endif()

message("CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")

message("Leaving geosxOptions.cmake\n")
