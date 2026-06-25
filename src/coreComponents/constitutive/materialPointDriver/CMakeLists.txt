# SPDX-License-Identifier: LGPL-2.1-only
#
# Optional, isolated MPM constitutive material-point driver support.
# Nothing in this directory is built unless GEOS_ENABLE_MPM_MATERIAL_POINT_DRIVER is ON.

set( geos_mpm_material_point_driver_headers
     MaterialPointDriver.hpp )

set( geos_mpm_material_point_driver_sources
     MaterialPointDriver.cpp )

blt_add_library( NAME       mpmMaterialPointDriver
                 SOURCES    ${geos_mpm_material_point_driver_sources}
                 HEADERS    ${geos_mpm_material_point_driver_headers}
                 DEPENDS_ON constitutive ${parallelDeps}
                 OBJECT     ${GEOS_BUILD_OBJ_LIBS}
                 SHARED     ${GEOS_BUILD_SHARED_LIBS} )

target_include_directories( mpmMaterialPointDriver PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )

blt_add_executable( NAME       geos_mpm_material_point_driver
                    SOURCES    MPMMaterialPointDriver.cpp
                    DEPENDS_ON constitutive ${parallelDeps} )

target_include_directories( geos_mpm_material_point_driver PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )

install( TARGETS mpmMaterialPointDriver LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib )
install( TARGETS geos_mpm_material_point_driver RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin )
