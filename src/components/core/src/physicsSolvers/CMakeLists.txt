
#
# Specify all headers
#
set( physicsSolvers_headers
     PhysicsSolverManager.hpp
     SolverBase.hpp
     SimpleSolvers/LaplaceFEM.hpp
     FiniteVolume/FlowSolverBase.hpp
     FiniteVolume/SinglePhaseFlow.hpp
     FiniteVolume/CompositionalMultiphaseFlow.hpp
     src/SolidMechanicsLagrangianFEM.hpp )

#
# Specify all sources
#
set( physicsSolvers_sources
     PhysicsSolverManager.cpp
     SolverBase.cpp
     SimpleSolvers/LaplaceFEM.cpp
     FiniteVolume/FlowSolverBase.cpp
     FiniteVolume/SinglePhaseFlow.cpp
     FiniteVolume/CompositionalMultiphaseFlow.cpp
     src/SolidMechanicsLagrangianFEM.cpp )


blt_add_library( NAME                  physicsSolvers
                 SOURCES               ${physicsSolvers_sources}
                 HEADERS               ${physicsSolvers_headers}
                 DEPENDS_ON            dataRepository systemSolverInterface constitutive finiteElement finiteVolume
               )
               
target_include_directories( physicsSolvers PUBLIC ${CMAKE_SOURCE_DIR}/components/core/src )
geosx_add_code_checks( PREFIX physicsSolvers )
