set( python_headers )

set( python_sources
     pygeosx.cpp )

set( dependencyList ${extraComponentsLinkList} ${externalComponentsLinkList} Python3::Python )

blt_add_library( NAME pygeosx
                 SOURCES ${python_sources}
                 HEADERS ${python_headers}
                 DEPENDS_ON ${dependencyList}
                 SHARED TRUE
                 CLEAR_PREFIX TRUE
                 )

target_include_directories( pygeosx PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )

geosx_add_code_checks( PREFIX pygeosx )
