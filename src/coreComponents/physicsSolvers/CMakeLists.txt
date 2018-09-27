
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
     src/SolidMechanicsLagrangianFEM.hpp
     surfaceGeneration/SurfaceGenerator.hpp
     src/DummySolver.hpp
     )

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
     src/SolidMechanicsLagrangianFEM.cpp
     surfaceGeneration/SurfaceGenerator.cpp
     src/DummySolver.cpp
     )


blt_add_library( NAME                  physicsSolvers
                 SOURCES               ${physicsSolvers_sources}
                 HEADERS               ${physicsSolvers_headers}
                 DEPENDS_ON            dataRepository systemSolverInterface
               )
               
target_include_directories( physicsSolvers PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )

geosx_add_code_checks( PREFIX physicsSolvers )

add_subdirectory( unitTests )