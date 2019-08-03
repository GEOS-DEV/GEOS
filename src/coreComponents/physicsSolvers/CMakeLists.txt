message(STATUS "Entering src/coreComponents/physicsSolvers/CMakeLists.txt")

#
# Specify solver headers
#
set( physicsSolvers_headers
     PhysicsSolverManager.hpp
     SolverBase.hpp
     SystemSolverParameters.hpp
     CoupledSolvers/HydrofractureSolver.hpp
     CoupledSolvers/PoroelasticSolver.hpp
     SimpleSolvers/LaplaceFEM.hpp
     FiniteVolume/FlowSolverBase.hpp
     FiniteVolume/SinglePhaseFlow.hpp
     FiniteVolume/SinglePhaseFlowKernels.hpp
     FiniteVolume/CompositionalMultiphaseFlow.hpp
     FiniteVolume/CompositionalMultiphaseFlowKernels.hpp
     solidMechanics/SolidMechanicsLagrangianFEM.hpp
     solidMechanics/SolidMechanicsLagrangianSSLE.hpp
     surfaceGeneration/SurfaceGenerator.hpp
     src/DummySolver.hpp
     solidMechanics/SolidMechanicsLagrangianFEMKernels.hpp
     solidMechanics/SolidMechanicsLagrangianSSLEKernels.hpp
     )

#
# Specify solver sources
#
set( physicsSolvers_sources
     PhysicsSolverManager.cpp
     SolverBase.cpp
     SystemSolverParameters.cpp
     CoupledSolvers/HydrofractureSolver.cpp
     CoupledSolvers/PoroelasticSolver.cpp
     SimpleSolvers/LaplaceFEM.cpp
     FiniteVolume/FlowSolverBase.cpp
     FiniteVolume/SinglePhaseFlow.cpp
     FiniteVolume/SinglePhaseFlowKernels.cpp
     FiniteVolume/CompositionalMultiphaseFlow.cpp
     FiniteVolume/CompositionalMultiphaseFlowKernels.cpp
     solidMechanics/SolidMechanicsLagrangianFEM.cpp
     solidMechanics/SolidMechanicsLagrangianSSLE.cpp
     surfaceGeneration/SurfaceGenerator.cpp
     src/DummySolver.cpp
     )

set( dependencyList common constitutive dataRepository linearAlgebraInterface )

set( externalComponentDeps "" )
if( ENABLE_GEOSX_PTP )
  message(STATUS "Adding GEOSX_PTP to physics solver dependencies")
  list( APPEND externalComponentDeps GEOSX_PTP )  
endif()

if ( ENABLE_OPENMP )
  set( dependencyList ${dependencyList} openmp )
endif()

if ( ENABLE_CUDA )
  set( dependencyList ${dependencyList} cuda )
endif()

blt_add_library( NAME                  physicsSolvers
                 SOURCES               ${physicsSolvers_sources}
                 HEADERS               ${physicsSolvers_headers}
                 DEPENDS_ON            ${dependencyList} ${externalComponentDeps}
               )

target_include_directories( physicsSolvers PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )
if( externalComponentDeps )
  target_include_directories( physicsSolvers PUBLIC ${CMAKE_SOURCE_DIR}/externalComponents)
endif()

geosx_add_code_checks( PREFIX physicsSolvers )

add_subdirectory( unitTests )
message(STATUS "Leaving src/coreComponents/physicsSolvers/CMakeLists.txt")
