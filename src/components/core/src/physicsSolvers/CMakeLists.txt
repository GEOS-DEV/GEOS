
#
# Specify all headers
#
set(physicsSolvers_headers
    PhysicsSolverManager.hpp
  	SolverBase.hpp
    BoundaryConditions/BoundaryConditionBase.hpp
    BoundaryConditions/BoundaryConditionManager.hpp
    BoundaryConditions/DirichletBoundaryCondition.hpp
    SimpleSolvers/LaplaceFEM.hpp )

#
# Specify all sources
#
set(physicsSolvers_sources
    PhysicsSolverManager.cpp
    SolverBase.cpp
    BoundaryConditions/BoundaryConditionBase.cpp
    BoundaryConditions/BoundaryConditionManager.cpp
    BoundaryConditions/DirichletBoundaryCondition.cpp
    SimpleSolvers/LaplaceFEM.cpp  )

#set( dependencyList dataRepository systemSolverInterface  )
#set( dependencyList2 "" )
#foreach( lib ${dependencyList} )
#    set( dependencyList2 ${dependencyList2} ${GEOSX_LINK_PREPEND_FLAG} ${lib} ${GEOSX_LINK_POSTPEND_FLAG} )
#endforeach()
#message( "dependencyList2 = ${dependencyList2}" )
#
#
#blt_add_library( NAME                  physicsSolvers
#                 SOURCES               ${physicsSolvers_sources}
#                 HEADERS               ${physicsSolvers_headers}
#                 DEPENDS_ON            ${dependencyList2}
#                 SHARED FALSE
#               )

blt_add_library( NAME                  physicsSolvers
                 SOURCES               ${physicsSolvers_sources}
                 HEADERS               ${physicsSolvers_headers}
                 DEPENDS_ON            dataRepository systemSolverInterface
                 SHARED FALSE
                 OBJECT TRUE
               )
                              
target_include_directories( physicsSolvers PUBLIC ${CMAKE_SOURCE_DIR}/components/core/src)
geosx_add_code_checks(PREFIX physicsSolvers )