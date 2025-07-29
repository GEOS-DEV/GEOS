##------------------------------------------------------------------------------
## add_third_party_libraries()
##
## Searches for the tpls.cmake file and includes it if found
## If this is not found, the assumption is that this is a compilation of the
## third party libraries.
##
##------------------------------------------------------------------------------
macro(add_third_party_libraries)
    find_file(TPLS_CMAKE tpls.cmake
        PATHS ${CMAKE_SOURCE_DIR}/../host-configs
              ${CMAKE_CURRENT_LIST_DIR}
              ${CMAKE_CURRENT_LIST_DIR}/..  )
    
    if ( TPLS_CMAKE )
        message( STATUS "tpls.cmake found: ${TPLS_CMAKE}")
        if ( NOT DEFINED GEOS_TPL_DIR )
            if ( DEFINED ENV{GEOS_TPL_DIR} )
                set( GEOS_TPL_DIR "$ENV{GEOS_TPL_DIR}" CACHE PATH "" FORCE )  
            else()
                message( FATAL_ERROR "GEOS_TPL_DIR not defined." )
            endif()
        endif()
        include( ${TPLS_CMAKE} )
    else()
        message( STATUS "tpls.cmake not found")
        message( STATUS "Assuming building TPLs" )
    endif()
endmacro(add_third_party_libraries)

##------------------------------------------------------------------------------
## set_compiler_options()
##
## Sets the default/common compiler options.
##
##------------------------------------------------------------------------------
macro(set_compiler_options)

    set( COMMON_FLAGS  ""              )
    set( RELEASE_FLAGS "-O3 -DNDEBUG"  )
    set( DEBUG_FLAGS   "-O0 -g"        )

    set( CMAKE_C_FLAGS               ${COMMON_FLAGS}  CACHE STRING "" )
    set( CMAKE_CXX_FLAGS             ${COMMON_FLAGS}  CACHE STRING "" )
    set( CMAKE_Fortran_FLAGS         ${COMMON_FLAGS}  CACHE STRING "" )
    set( CMAKE_CXX_FLAGS_RELEASE     ${RELEASE_FLAGS} CACHE STRING "" )
    set( CMAKE_C_FLAGS_RELEASE       ${RELEASE_FLAGS} CACHE STRING "" )
    set( CMAKE_Fortran_FLAGS_RELEASE ${RELEASE_FLAGS} CACHE STRING "" )
    set( CMAKE_CXX_FLAGS_DEBUG       ${DEBUG_FLAGS}   CACHE STRING "" )
    set( CMAKE_C_FLAGS_DEBUG         ${DEBUG_FLAGS}   CACHE STRING "" )
    set( CMAKE_Fortran_FLAGS_DEBUG   ${DEBUG_FLAGS}   CACHE STRING "" )

endmacro(set_compiler_options)

##------------------------------------------------------------------------------
## set_cuda_options()
##
## Sets the default/common CUDA options.
##
##------------------------------------------------------------------------------
macro(set_cuda_options)

    set( CUDA_TOOLKIT_ROOT_DIR "$ENV{CUDA_HOME}" CACHE PATH "" )
    set( CMAKE_CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER} CACHE STRING "" )
    set( CMAKE_CUDA_COMPILER ${CUDA_TOOLKIT_ROOT_DIR}/bin/nvcc CACHE STRING "" )
    set( CUDA_ARCH "sm_90" CACHE STRING "" )
    set( CMAKE_CUDA_STANDARD 17 CACHE STRING "" )
    set( CMAKE_CUDA_FLAGS "-restrict -arch ${CUDA_ARCH} --expt-extended-lambda --expt-relaxed-constexpr -Werror cross-execution-space-call,reorder,deprecated-declarations -Xcompiler -std=c++17" CACHE STRING "" )
    set( CMAKE_CUDA_FLAGS_RELEASE "-O3 -DNDEBUG -Xcompiler -DNDEBUG -Xcompiler -O3" CACHE STRING "" )
    set( CMAKE_CUDA_FLAGS_RELWITHDEBINFO "-g -lineinfo ${CMAKE_CUDA_FLAGS_RELEASE}" CACHE STRING "" )
    set( CMAKE_CUDA_FLAGS_DEBUG "-g -G -O0 -Xcompiler -O0" CACHE STRING "" )

endmacro(set_cuda_options)

#########################################################################
# BLAS/LAPACK SETUP                                                                                                                                                               
#########################################################################
# use :
# - OpenBLAS library
find_library( BLAS_LIBRARIES 
    NAMES openblas openblas.so.0 libopenblas.so.0
    PATHS /usr
    PATH_SUFFIXES lib lib64 )

if(NOT BLAS_LIBRARIES)
    message( FATAL_ERROR "Failed to find OpenBLAS." )
endif()

set( LAPACK_LIBRARIES ${BLAS_LIBRARIES} )

message( STATUS "OpenBLAS BLAS '${BLAS_LIBRARIES}'." )
message( STATUS "OpenBLAS LAPACK '${LAPACK_LIBRARIES}'." )

set( ENABLE_MKL OFF CACHE BOOL "" )

#########################################################################
# OTHER OPTIONS                                                                                                                                                            
#########################################################################
# Use old policy CMP0146
set ( CMAKE_POLICY_DEFAULT_CMP0146 "OLD" CACHE STRING "" FORCE )

# Use new policy CMP0135
set ( CMAKE_POLICY_DEFAULT_CMP0135 "NEW" CACHE STRING "" FORCE )

# Enable failure on tests
set( ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "" FORCE )

# Enale caliper instrumentation
set( ENABLE_CALIPER ON CACHE BOOL "" )

# Disable doxygen. No valid libiconv available
set( ENABLE_DOXYGEN OFF CACHE BOOL "" FORCE )

# Disable MathPresso. Not supported in this architecture
set( ENABLE_MATHPRESSO OFF CACHE BOOL "" FORCE )

