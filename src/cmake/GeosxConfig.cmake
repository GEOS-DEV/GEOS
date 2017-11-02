#
# Get GEOSX Version
file (STRINGS "VERSION" GEOSX_VERSION_FULL)
string(REGEX REPLACE "VERSION_ID = v" "" GEOSX_VERSION_FULL "${GEOSX_VERSION_FULL}")
string(REPLACE "." ";" GEOSX_VERSION_LIST ${GEOSX_VERSION_FULL})

list( GET GEOSX_VERSION_LIST  0 GEOSX_VERSION_MAJOR )
list( GET GEOSX_VERSION_LIST  1 GEOSX_VERSION_MINOR )
list( GET GEOSX_VERSION_LIST  2 GEOSX_VERSION_PATCH )
message(STATUS "Configuring GEOSX version ${GEOSX_VERSION_FULL}")


set(TPL_DEPENDENCIES ATK CALIPER CHAI FPARSER MATHPRESSO PYTHON RAJA )

foreach( DEP in ${TPL_DEPENDENCIES})
    if( ${DEP}_FOUND OR ENABLE_${DEP} )
        set(USE_${DEP} TRUE  )
    endif()
endforeach()


configure_file(
    ${CMAKE_SOURCE_DIR}/components/core/src/common/GeosxConfig.hpp.in
    ${CMAKE_CURRENT_BINARY_DIR}/include/common/GeosxConfig.hpp
)
configure_file(
    ${CMAKE_SOURCE_DIR}/components/core/src/common/GeosxConfig.hpp.in
    ${CMAKE_SOURCE_DIR}/components/core/src/common/GeosxConfig.hpp
)

install( FILES ${CMAKE_CURRENT_BINARY_DIR}/include/common/GeosxConfig.hpp
         DESTINATION ${CMAKE_INSTALL_PREFIX}/include/common )

