#
# Specify all headers
#
set( mainInterface_headers
     GeosxState.hpp
     ProblemManager.hpp
     initialization.hpp
     version.hpp
     ControlledInput.hpp
   )

#
# Specify all sources
#
set( mainInterface_sources
     GeosxState.cpp
     ProblemManager.cpp
     initialization.cpp
     version.cpp
     ControlledInput.cpp
   )

set( dependencyList ${parallelDeps} physicsSolvers discretizationMethods fieldSpecification linearAlgebra dataRepository events fileIO optionparser )

if( ENABLE_CONTROLLED_INPUT )
  set( dependencyList ${dependencyList} yaml-cpp )
endif()

blt_add_library( NAME       mainInterface
                 SOURCES    ${mainInterface_sources}
                 HEADERS    ${mainInterface_headers}
                 DEPENDS_ON ${dependencyList}
                 OBJECT     ${GEOSX_BUILD_OBJ_LIBS}
               )

add_dependencies( mainInterface generate_version )

target_include_directories( mainInterface PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )

geosx_add_code_checks( PREFIX mainInterface )
