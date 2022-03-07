set( pyfileIO_headers PyHistoryCollectionType.hpp
  PyHistoryOutputType.hpp
		    )

set( pyfileIO_sources PyHistoryCollection.cpp
PyHistoryOutput.cpp
  )


blt_add_library( NAME             pyfileIO
                 SOURCES          ${pyfileIO_sources}
                 HEADERS          ${pyfileIO_headers}
                 DEPENDS_ON       ${dependencyList} ${extraComponentsList} ${externalComponentList} Python3::Python fileIO
                 OBJECT           ${GEOSX_BUILD_OBJ_LIBS}
                 )

target_include_directories( pyfileIO PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )

geosx_add_code_checks( PREFIX pyfileIO )
