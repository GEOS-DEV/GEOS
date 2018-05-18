
#
# Specify all headers
#
set( codingUtilities_headers
     StringUtilities.hpp
     IOUtilities.hpp
   )

#
# Specify all sources
#
set( codingUtilities_sources
     StringUtilities.cpp
   )



set( dependencyList common RAJA )

if( CONDUIT_FOUND )
  set( dependencyList ${dependencyList} conduit )
endif()

blt_add_library( NAME                  codingUtilities
                 SOURCES               ${codingUtilities_sources}
                 HEADERS               ${codingUtilities_headers}
                 DEPENDS_ON            ${dependencyList}
#                 SHARED TRUE
               )
               
target_include_directories( codingUtilities PUBLIC ${CMAKE_SOURCE_DIR}/components/core/src)
geosx_add_code_checks(PREFIX codingUtilities )