set( pyfileIO_headers PyHistoryCollectionType.hpp
  PyHistoryOutputType.hpp
		    )

set( pyfileIO_sources PyHistoryCollection.cpp
PyHistoryOutput.cpp
  )

set( dependencyList ${extraComponentsLinkList} ${externalComponentsLinkList} Python3::Python )

blt_add_library( NAME             pyfileIO
                 SOURCES          ${pyfileIO_sources}
                 HEADERS          ${pyfileIO_headers}
                 DEPENDS_ON       fileIO mainInterface ${dependencyList}
                 SHARED TRUE
                 CLEAR_PREFIX TRUE
                 )

target_include_directories( pyfileIO PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )

geosx_add_code_checks( PREFIX pyfileIO )
