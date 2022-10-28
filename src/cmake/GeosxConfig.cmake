set( PREPROCESSOR_DEFINES ARRAY_BOUNDS_CHECK
                          CALIPER
                          CHAI
                          CUDA
                          CUDA_NVTOOLSEXT
                          FORTRAN_MANGLE_NO_UNDERSCORE
                          FPE
                          HYPRE
                          HYPRE_CUDA
                          MATHPRESSO
                          METIS
                          MKL
                          MPI
                          OPENMP
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

    configure_file( ${CMAKE_SOURCE_DIR}/coreComponents/common/GeosxConfig.hpp.in
                    ${CMAKE_SOURCE_DIR}/docs/doxygen/GeosxConfig.hpp )
endfunction()


make_full_config_file( "${PREPROCESSOR_DEFINES}" )
