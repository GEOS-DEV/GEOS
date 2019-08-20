
#
# Specify all headers
#
set(systemSolverInterface_headers
    EpetraBlockSystem.hpp
    LinearSolverWrapper.hpp
    )

#
# Specify all sources
#
set(systemSolverInterface_sources
    EpetraBlockSystem.cpp
    LinearSolverWrapper.cpp
   )

set( dependencyList ${dependencyList} fileIO common trilinos )

blt_add_library( NAME                  systemSolverInterface
                 SOURCES               ${systemSolverInterface_sources}
                 HEADERS               ${systemSolverInterface_headers}
                 DEPENDS_ON            ${dependencyList}
               )
               
target_include_directories( systemSolverInterface PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)

geosx_add_code_checks(PREFIX systemSolverInterface )
