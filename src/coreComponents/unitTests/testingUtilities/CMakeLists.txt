#
# Specify all headers
#
set( testingUtilities_headers
     TestingTasks.hpp
   )

#
# Specify all sources
#
set( testingUtilities_sources
     TestingTasks.cpp
   )

#
# Specify all dependencies
#
set( dependencyList ${parallelDeps} gtest )

if ( GEOS_BUILD_SHARED_CORE_LIB )
  list( APPEND dependencyList geosx_core HDF5::HDF5 )
else()
  list( APPEND dependencyList ${geosx_core_libs} HDF5::HDF5 )
endif()
   
blt_add_library( NAME       testingUtilities
                 SOURCES    ${testingUtilities_sources}
                 HEADERS    ${testingUtilities_headers}
                 DEPENDS_ON ${dependencyList}
                 OBJECT     ${GEOS_BUILD_OBJ_LIBS}
                 SHARED     ${GEOS_BUILD_SHARED_LIBS}
               )

target_include_directories( testingUtilities PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )
