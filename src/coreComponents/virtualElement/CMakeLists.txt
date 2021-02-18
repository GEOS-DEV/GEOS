#
# Specify all headers
#
set( virtualElement_headers
  ConformingVirtualElementOrder1.hpp
  VirtualElementBase.hpp
  )
#
# Specify all sources
#
set( virtualElement_sources
  ConformingVirtualElementOrder1.cpp
  )

if( BUILD_OBJ_LIBS)
  set( dependencyList dataRepository codingUtilities)
else()
  set( dependencyList mesh)
endif()

if( ENABLE_OPENMP )
    set( dependencyList ${dependencyList} openmp )
endif()

if ( ENABLE_CUDA )
  set( dependencyList ${dependencyList} cuda )
endif()

blt_add_library( NAME                  virtualElement
                 SOURCES               ${virtualElement_sources}
                 HEADERS               ${virtualElement_headers}
                 DEPENDS_ON            ${dependencyList}
                 OBJECT                ${buildAsObj}
               )

target_include_directories( virtualElement PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)

geosx_add_code_checks( PREFIX virtualElement )

add_subdirectory( unitTests )
