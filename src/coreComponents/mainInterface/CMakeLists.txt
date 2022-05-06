#
# Specify all headers
#
set( mainInterface_headers
     GeosxState.hpp
     ProblemManager.hpp
     initialization.hpp
     StableToXml.hpp
   )

#
# Specify all sources
#
set( mainInterface_sources
     GeosxState.cpp
     ProblemManager.cpp
     initialization.cpp
     StableToXml.cpp
   )

set( dependencyList physicsSolvers discretizationMethods fieldSpecification linearAlgebra events fileIO optionparser )

if( ENABLE_STABILIZED_INPUT )
  set( dependencyList ${dependencyList} yaml-cpp )
endif()

if ( ENABLE_CUDA )
  set( dependencyList ${dependencyList} cuda )
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
