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
                 DEPENDS_ON       ${dependencyList} ${extraComponentsList} ${externalComponentList} Python3::Python dataRepository mainInterface
		 OBJECT           ${GEOSX_BUILD_OBJ_LIBS}                 )

target_include_directories( pygroup PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )

geosx_add_code_checks( PREFIX pygroup )
