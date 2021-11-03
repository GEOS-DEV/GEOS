set( python_headers pygeosx.hpp
                    PyGroup.hpp
                    PyWrapper.hpp
                    PySolver.hpp
		    PyWaveSolver.hpp
                    )

set( python_sources pygeosx.cpp
                    PyGroup.cpp
                    PyWrapper.cpp
                    PySolver.cpp
		    PyWaveSolver.cpp
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
