
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

set( dependencyList common )
set( dependencyList2 RAJA )
foreach( lib ${dependencyList} )
    set( dependencyList2 ${dependencyList2} ${GEOSX_LINK_PREPEND_FLAG} ${lib} ${GEOSX_LINK_POSTPEND_FLAG} )
endforeach()
message( "dependencyList2 = ${dependencyList2}" )
    
blt_add_library( NAME                  codingUtilities
                 SOURCES               ${codingUtilities_sources}
                 HEADERS               ${codingUtilities_headers}
                 DEPENDS_ON            ${dependencyList2}
                 SHARED FALSE
               )
               
               
               
target_include_directories( codingUtilities PUBLIC ${CMAKE_SOURCE_DIR}/components/core/src)
geosx_add_code_checks(PREFIX codingUtilities )