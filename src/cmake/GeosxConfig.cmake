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
                          CONTAINERARRAY_RETURN_PTR
                          FPARSER
                          HYPRE
                          MATHPRESSO
                          METIS
                          MPI
                          OPENMP
                          PARMETIS
                          PYTHON
                          RAJA 
                          SUPERLU_DIST
                          TIMERS
                          TRILINOS
                          ${externalComponentsList} )

foreach( DEP in ${PREPROCESSOR_DEFINES})
    if( ${DEP}_FOUND OR ENABLE_${DEP} )
        set(USE_${DEP} TRUE  )
        set(GEOSX_USE_${DEP} TRUE  )
    endif()
endforeach()


configure_file(
    ${CMAKE_SOURCE_DIR}/coreComponents/common/GeosxConfig.hpp.in
    ${CMAKE_BINARY_DIR}/include/common/GeosxConfig.hpp
)


# This approach requires a complete rebuild when switching between different build configurations,
# as a new GeosxConfig.hpp file is created when you switch builds.
#configure_file(
#    ${CMAKE_SOURCE_DIR}/coreComponents/common/GeosxConfig.hpp.in
#    ${CMAKE_SOURCE_DIR}/coreComponents/common/GeosxConfig.hpp
#)

# This approach does not. I guess the symbolic link points to the file in the build directory, and
# uses that date. So it only rebuilds files that have changed since the last build in the specific
# configuration.
#ADD_CUSTOM_TARGET(geosx_config_hpp ALL
#                  COMMAND ${CMAKE_COMMAND} -E create_symlink
#                  ${CMAKE_BINARY_DIR}/include/common/GeosxConfig.hpp
#                  ${CMAKE_SOURCE_DIR}/coreComponents/common/GeosxConfig.hpp )


install( FILES ${CMAKE_BINARY_DIR}/include/common/GeosxConfig.hpp
         DESTINATION ${CMAKE_INSTALL_PREFIX}/include/common )

