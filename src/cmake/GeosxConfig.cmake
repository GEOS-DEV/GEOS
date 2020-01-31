#
# Get GEOSX Version
file (STRINGS "VERSION" GEOSX_VERSION_FULL)
string(REGEX REPLACE "VERSION_ID = v" "" GEOSX_VERSION_FULL "${GEOSX_VERSION_FULL}")
string(REPLACE "." ";" GEOSX_VERSION_LIST ${GEOSX_VERSION_FULL})

list( GET GEOSX_VERSION_LIST  0 GEOSX_VERSION_MAJOR )
list( GET GEOSX_VERSION_LIST  1 GEOSX_VERSION_MINOR )
list( GET GEOSX_VERSION_LIST  2 GEOSX_VERSION_PATCH )
message(STATUS "Configuring GEOSX version ${GEOSX_VERSION_FULL}")


set( PREPROCESSOR_DEFINES ARRAY_BOUNDS_CHECK
                          ATK
                          CALIPER
                          CHAI
                          CUDA
                          FORTRAN_MANGLE_NO_UNDERSCORE
                          HYPRE
                          MATHPRESSO
                          METIS
                          MPI
                          OPENMP
                          CUDA
                          PARMETIS
                          PETSC
                          PYTHON
                          RAJA 
                          SUPERLU_DIST
                          TIMERS
                          TOTALVIEW_OUTPUT
                          TRILINOS
                          MKL
                          GEOSX_PTP
                          SEPARATION_COEFFICIENT
                          ${externalComponentsList} )

foreach( DEP in ${PREPROCESSOR_DEFINES})
    if( ${DEP}_FOUND OR ENABLE_${DEP} )
        set(USE_${DEP} TRUE  )
        set(GEOSX_USE_${DEP} TRUE  )
    endif()
endforeach()


configure_file( ${CMAKE_SOURCE_DIR}/coreComponents/common/GeosxConfig.hpp.in
                ${CMAKE_BINARY_DIR}/include/common/GeosxConfig.hpp )

install( FILES ${CMAKE_BINARY_DIR}/include/common/GeosxConfig.hpp
         DESTINATION ${CMAKE_INSTALL_PREFIX}/include/common )


function( make_full_config_file 
          PREPROCESSOR_VARS )
    foreach( DEP in ${PREPROCESSOR_VARS})
        set(USE_${DEP} TRUE  )
        set(GEOSX_USE_${DEP} TRUE  )
        set(${DEP} TRUE  )
    endforeach()

    configure_file( ${CMAKE_SOURCE_DIR}/coreComponents/common/GeosxConfig.hpp.in
                    ${CMAKE_SOURCE_DIR}/docs/doxygen/GeosxConfig.hpp )
endfunction()


make_full_config_file( "${PREPROCESSOR_DEFINES}" )
