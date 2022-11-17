set( PREPROCESSOR_DEFINES ARRAY_BOUNDS_CHECK
                          CALIPER
                          CHAI
                          CUDA
                          CUDA_NVTOOLSEXT
                          HIP
                          FORTRAN_MANGLE_NO_UNDERSCORE
                          FPE
                          HYPRE
                          MATHPRESSO
                          METIS
                          MKL
                          MPI
                          PARMETIS
                          PETSC
                          PVTPackage
                          PYGEOSX
                          RAJA
                          SCOTCH
                          SEPARATION_COEFFICIENT
                          SUITESPARSE
                          SUPERLU_DIST
                          TIMERS
                          TOTALVIEW_OUTPUT
                          TRILINOS
                          VTK
                          ${externalComponentsList} )

foreach( DEP in ${PREPROCESSOR_DEFINES} )
    if( ${DEP}_FOUND OR ENABLE_${DEP} OR GEOSX_ENABLE_${DEP} )
        set( USE_${DEP} TRUE )
        set( GEOSX_USE_${DEP} TRUE )
	message(STATUS "GEOSX_USE_${DEP} = ${GEOSX_USE_${DEP}}")
    endif()
endforeach()

set( STRICT_PPD OPENMP )

# only activate these options if they are ENABLED AND FOUND, not if either
foreach( DEP in ${STRICT_PPD} )
    if( ${DEP}_FOUND AND ( ENABLE_${DEP} OR GEOSX_ENABLE_${DEP} ) )
        set( USE_${DEP} TRUE )
        set( GEOSX_USE_${DEP} TRUE )
	message(STATUS "GEOSX_USE_${DEP} = ${GEOSX_USE_${DEP}}")
    endif()
endforeach( )

set( HYPRE_DEVICE_OPTIONS CPU CUDA HIP )
if( NOT ENABLE_HYPRE_DEVICE )
  set( ENABLE_HYPRE_DEVICE CPU )
endif()
if(NOT ${ENABLE_HYPRE_DEVICE} IN_LIST HYPRE_DEVICE_OPTIONS )
    message(FATAL_ERROR "Set ENABLE_HYPRE_DEVICE to CPU, CUDA, or HIP.")
endif()
set( GEOSX_USE_HYPRE_DEVICE "GEOSX_USE_HYPRE_${ENABLE_HYPRE_DEVICE}" )
message( STATUS "GEOSX_USE_HYPRE_DEVICE = ${GEOSX_USER_HYPRE_DEVICE}")

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

    configure_file( ${CMAKE_SOURCE_DIR}/coreComponents/common/GeosxConfig.hpp.in
                    ${CMAKE_SOURCE_DIR}/docs/doxygen/GeosxConfig.hpp )
endfunction()


make_full_config_file( "${PREPROCESSOR_DEFINES}" )
