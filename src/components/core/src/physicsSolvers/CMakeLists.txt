
#
# Specify all headers
#
set(physicsSolvers_headers
    PhysicsSolverManager.hpp
    SolverBase.hpp
    BoundaryConditions/BoundaryConditionBase.hpp
    BoundaryConditions/BoundaryConditionManager.hpp
    BoundaryConditions/DirichletBoundaryCondition.hpp
    FiniteVolume/SinglePhaseFlow_TPFA.hpp
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
    FiniteVolume/SinglePhaseFlow_TPFA.cpp
    SimpleSolvers/LaplaceFEM.cpp  )


blt_add_library( NAME                  physicsSolvers
                 SOURCES               ${physicsSolvers_sources}
                 HEADERS               ${physicsSolvers_headers}
                 DEPENDS_ON            dataRepository systemSolverInterface
               )
               
target_include_directories( physicsSolvers PUBLIC ${CMAKE_SOURCE_DIR}/components/core/src)
geosx_add_code_checks(PREFIX physicsSolvers )