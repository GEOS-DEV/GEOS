set( pygroup_headers PyGroup.hpp
                     PyGroupType.hpp
                     PyWrapper.hpp
		    )

set( pygroup_sources PyGroup.cpp
                     PyWrapper.cpp
  )


blt_add_library( NAME             pygroup
                 SOURCES          ${pygroup_sources}
                 HEADERS          ${pygroup_headers}
                 DEPENDS_ON       ${dependencyList} ${extraComponentsLinkList} ${externalComponentsLinkList} Python3::Python dataRepository mainInterface
                 SHARED TRUE
                 CLEAR_PREFIX TRUE
                 )

target_include_directories( pygroup PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )

geosx_add_code_checks( PREFIX pygroup )
