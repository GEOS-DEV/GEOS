
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


#set( dependencyList fileIO common )
#set( dependencyList2 trilinos )
#foreach( lib ${dependencyList} )
#    set( dependencyList2 ${dependencyList2} ${GEOSX_LINK_PREPEND_FLAG} ${lib} ${GEOSX_LINK_POSTPEND_FLAG} )
#endforeach()
#message( "dependencyList2 = ${dependencyList2}" )

blt_add_library( NAME                  systemSolverInterface
                 SOURCES               ${systemSolverInterface_sources}
                 HEADERS               ${systemSolverInterface_headers}
                 DEPENDS_ON            fileIO common trilinos
                 SHARED FALSE
               )
               
target_include_directories( systemSolverInterface PUBLIC ${CMAKE_SOURCE_DIR}/components/core/src)
geosx_add_code_checks(PREFIX systemSolverInterface )