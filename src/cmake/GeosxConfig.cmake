#
# Get GEOSX Version
#
file ( STRINGS "VERSION" GEOSX_VERSION_FULL )
string( REGEX REPLACE "VERSION_ID = v" "" GEOSX_VERSION_FULL "${GEOSX_VERSION_FULL}" )
string( REPLACE "." ";" GEOSX_VERSION_LIST ${GEOSX_VERSION_FULL} )

list( GET GEOSX_VERSION_LIST  0 GEOSX_VERSION_MAJOR )
list( GET GEOSX_VERSION_LIST  1 GEOSX_VERSION_MINOR )
list( GET GEOSX_VERSION_LIST  2 GEOSX_VERSION_PATCH )

if( GIT_FOUND )
    blt_is_git_repo( OUTPUT_STATE is_git_repo )
    if( is_git_repo )
        # the other option is to use blt_git_tag() but we want branch info too
        blt_git_branch( BRANCH_NAME git_branch RETURN_CODE git_branch_rc )
        blt_git_hashcode( HASHCODE git_hash RETURN_CODE git_hash_rc )
        set( GEOSX_VERSION_DEV "${git_branch}-${git_hash}" )
    endif()
endif()

message( STATUS "Configuring GEOSX version ${GEOSX_VERSION_FULL} ${GEOSX_VERSION_DEV}" )

set( PREPROCESSOR_DEFINES ARRAY_BOUNDS_CHECK
                          CALIPER
                          CHAI
                          CUDA
                          FORTRAN_MANGLE_NO_UNDERSCORE
                          FPE
                          HYPRE
                          HYPRE_CUDA
                          MATHPRESSO
                          METIS
                          MPI
                          OPENMP
                          CUDA
                          PARMETIS
                          PETSC
                          PVTPackage
                          PYGEOSX
                          RAJA 
                          SUPERLU_DIST
                          SUITESPARSE
                          TIMERS
                          TOTALVIEW_OUTPUT
                          TRILINOS
                          MKL
                          SEPARATION_COEFFICIENT
                          ${externalComponentsList} )

foreach( DEP in ${PREPROCESSOR_DEFINES} )
    if( ${DEP}_FOUND OR ENABLE_${DEP} OR GEOSX_ENABLE_${DEP} )
        set( USE_${DEP} TRUE )
        set( GEOSX_USE_${DEP} TRUE )
    endif()
endforeach()

set( GEOSX_CMAKE_BUILD_TYPE "\"${CMAKE_BUILD_TYPE}\"" )

configure_file( ${CMAKE_SOURCE_DIR}/coreComponents/common/GeosxConfig.hpp.in
                ${CMAKE_BINARY_DIR}/include/common/GeosxConfig.hpp )

install( FILES ${CMAKE_BINARY_DIR}/include/common/GeosxConfig.hpp
         DESTINATION ${CMAKE_INSTALL_PREFIX}/include/common )


function( make_full_config_file 
          PREPROCESSOR_VARS )
    foreach( DEP in ${PREPROCESSOR_VARS} )
        set( USE_${DEP} TRUE )
        set( GEOSX_USE_${DEP} TRUE )
        set( ${DEP} TRUE )
    endforeach()

    # Fix some options to avoid changes in committed config for doxygen
    set( GEOSX_CMAKE_BUILD_TYPE "\"Release\"" )
    set( GEOSX_LOCALINDEX_TYPE "std::ptrdiff_t" )
    set( GEOSX_LOCALINDEX_TYPE_FLAG "3" )
    set( GEOSX_GLOBALINDEX_TYPE "long long int" )
    set( GEOSX_GLOBALINDEX_TYPE_FLAG "2" )
    set( GEOSX_LA_INTERFACE "Hypre" )
    set( GEOSX_LA_INTERFACE_HYPRE ON )
    set( GEOSX_LA_INTERFACE_TRILINOS OFF )
    set( GEOSX_LA_INTERFACE_PETSC OFF )
    unset( GEOSX_VERSION_DEV )

    configure_file( ${CMAKE_SOURCE_DIR}/coreComponents/common/GeosxConfig.hpp.in
                    ${CMAKE_SOURCE_DIR}/docs/doxygen/GeosxConfig.hpp )
endfunction()


make_full_config_file( "${PREPROCESSOR_DEFINES}" )
