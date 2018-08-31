
#
# Specify all headers
#
set(systemSolverInterface_headers
    EpetraBlockSystem.hpp
    LinearSolverWrapper.hpp
    SystemSolverParameters.hpp
    )

#
# Specify all sources
#
set(systemSolverInterface_sources
    EpetraBlockSystem.cpp
    SystemSolverParameters.cpp
    LinearSolverWrapper.cpp
   )

blt_add_library( NAME                  systemSolverInterface
                 SOURCES               ${systemSolverInterface_sources}
                 HEADERS               ${systemSolverInterface_headers}
                 DEPENDS_ON            fileIO common trilinos
               )
               
target_include_directories( systemSolverInterface PUBLIC ${CMAKE_SOURCE_DIR}/components/core/src)
target_include_directories( systemSolverInterface PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)

geosx_add_code_checks(PREFIX systemSolverInterface )