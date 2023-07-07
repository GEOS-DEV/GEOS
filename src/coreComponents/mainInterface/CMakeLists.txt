#
# Specify all headers
#
set( mainInterface_headers
     GeosxState.hpp
     ProblemManager.hpp
     initialization.hpp
     version.hpp
     reproducer.hpp
   )

#
# Specify all sources
#
set( mainInterface_sources
    #  GeosxState.cpp
    #  ProblemManager.cpp
    #  initialization.cpp
    #  version.cpp
     reproducer.cpp
   )

set( dependencyList ${parallelDeps} dataRepository optionparser )

blt_add_library( NAME       mainInterface
                 SOURCES    ${mainInterface_sources}
                 HEADERS    ${mainInterface_headers}
                 DEPENDS_ON ${dependencyList}
                 OBJECT     ${GEOSX_BUILD_OBJ_LIBS}
               )

add_dependencies( mainInterface generate_version )

target_include_directories( mainInterface PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )

geosx_add_code_checks( PREFIX mainInterface )
