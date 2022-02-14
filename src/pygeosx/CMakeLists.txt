
set( python_headers pygeosx.hpp
                    ${CMAKE_SOURCE_DIR}/coreComponents/dataRepository/python/PyGroup.hpp
                    ${CMAKE_SOURCE_DIR}/coreComponents/dataRepository/python/PyGroupType.hpp
                    ${CMAKE_SOURCE_DIR}/coreComponents/dataRepository/python/PyWrapper.hpp
		    ${CMAKE_SOURCE_DIR}/coreComponents/fileIO/python/PyHistoryCollectionType.hpp
		    ${CMAKE_SOURCE_DIR}/coreComponents/fileIO/python/PyHistoryOutputType.hpp
		    ${CMAKE_SOURCE_DIR}/coreComponents/physicsSolvers/python/PySolverType.hpp
		    ${CMAKE_SOURCE_DIR}/coreComponents/physicsSolvers/python/PySolver.hpp
		    )

set( python_sources pygeosx.cpp
		    ${CMAKE_SOURCE_DIR}/coreComponents/dataRepository/python/PyGroup.cpp
                    ${CMAKE_SOURCE_DIR}/coreComponents/dataRepository/python/PyWrapper.cpp
		    ${CMAKE_SOURCE_DIR}/coreComponents/fileIO/python/PyHistoryCollection.cpp
		    ${CMAKE_SOURCE_DIR}/coreComponents/fileIO/python/PyHistoryOutput.cpp
		    ${CMAKE_SOURCE_DIR}/coreComponents/physicsSolvers/python/PySolver.cpp
		    )

set( dependencyList ${extraComponentsLinkList} ${externalComponentsLinkList} Python3::Python )

blt_add_library( NAME pygeosx
                 SOURCES ${python_sources}
                 HEADERS ${python_headers}
                 DEPENDS_ON ${dependencyList}
                 SHARED TRUE
                 CLEAR_PREFIX TRUE
               )

add_custom_command( TARGET pygeosx
                    PRE_BUILD
                    COMMAND ${Python3_EXECUTABLE} -m venv ${CMAKE_BINARY_DIR}/lib/PYGEOSX --system-site-packages )

add_custom_command( TARGET pygeosx
                    POST_BUILD
                    COMMAND ln -sf ${CMAKE_BINARY_DIR}/lib/pygeosx.so ${CMAKE_BINARY_DIR}/lib/PYGEOSX/lib/python3.8/site-packages/ )

add_custom_command( TARGET pygeosx
                    POST_BUILD
                    COMMAND ln -sf ${CMAKE_BINARY_DIR}/lib/pylvarray.so ${CMAKE_BINARY_DIR}/lib/PYGEOSX/lib/python3.8/site-packages/ )

target_include_directories( pygeosx PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )

geosx_add_code_checks( PREFIX pygeosx )
