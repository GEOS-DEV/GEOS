message( "Entering src/coreComponents/physicsSolvers/CMakeLists.txt")

#
# Specify all headers
#
set( physicsSolvers_headers
     PhysicsSolverManager.hpp
     SolverBase.hpp
     CoupledSolvers/PoroelasticSolver.hpp
     CoupledSolvers/ReservoirWellsSystemSolver.hpp
     SimpleSolvers/LaplaceFEM.hpp
     FiniteVolume/FlowSolverBase.hpp
     FiniteVolume/SinglePhaseFlow.hpp
     FiniteVolume/CompositionalMultiphaseFlow.hpp
     Wells/WellSolverBase.hpp
     Wells/SinglePhaseWell.hpp
     Wells/CompositionalMultiphaseWell.hpp
     solidMechanics/SolidMechanicsLagrangianFEM.hpp
     solidMechanics/SolidMechanicsLagrangianSSLE.hpp
     surfaceGeneration/SurfaceGenerator.hpp
     src/DummySolver.hpp
     )

#
# Specify all sources
#
set( physicsSolvers_sources
     PhysicsSolverManager.cpp
     SolverBase.cpp
     CoupledSolvers/PoroelasticSolver.cpp
     CoupledSolvers/ReservoirWellsSystemSolver.cpp
     SimpleSolvers/LaplaceFEM.cpp
     FiniteVolume/FlowSolverBase.cpp
     FiniteVolume/SinglePhaseFlow.cpp
     FiniteVolume/CompositionalMultiphaseFlow.cpp
     Wells/WellSolverBase.cpp
     Wells/SinglePhaseWell.cpp
     Wells/CompositionalMultiphaseWell.cpp
     solidMechanics/SolidMechanicsLagrangianFEM.cpp
     solidMechanics/SolidMechanicsLagrangianSSLE.cpp
     surfaceGeneration/SurfaceGenerator.cpp
     src/DummySolver.cpp
     )

set( externalComponentDeps "" )
if( ENABLE_GEOSX_PTP )
  message("adding GEOSX_PTP to physics solver dependencies")
  list( APPEND externalComponentDeps GEOSX_PTP )  
endif()

blt_add_library( NAME                  physicsSolvers
                 SOURCES               ${physicsSolvers_sources}
                 HEADERS               ${physicsSolvers_headers}
                 DEPENDS_ON            dataRepository systemSolverInterface ${externalComponentDeps}
               )
               
target_include_directories( physicsSolvers PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )
if( externalComponentDeps )
  target_include_directories( physicsSolvers PUBLIC ${CMAKE_SOURCE_DIR}/externalComponents)
endif()

geosx_add_code_checks( PREFIX physicsSolvers )

add_subdirectory( unitTests )
message( "Leaving src/coreComponents/physicsSolvers/CMakeLists.txt")
