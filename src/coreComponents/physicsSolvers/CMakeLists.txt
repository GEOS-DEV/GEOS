message(STATUS "Entering src/coreComponents/physicsSolvers/CMakeLists.txt")

add_subdirectory (GEOSX_PTP)

#
# Specify solver headers
#
set( physicsSolvers_headers
     PhysicsSolverManager.hpp
     SolverBase.hpp
     SystemSolverParameters.hpp
     NonlinearSolverParameters.hpp
     fluidFlow/CompositionalMultiphaseFlow.hpp
     fluidFlow/CompositionalMultiphaseFlowKernels.hpp
     fluidFlow/CompositionalMultiphaseWell.hpp
     fluidFlow/FlowSolverBase.hpp
     fluidFlow/SinglePhaseBase.hpp
     fluidFlow/SinglePhaseHybridFVM.hpp
     fluidFlow/SinglePhaseFVM.hpp
     fluidFlow/SinglePhaseKernels.hpp
     fluidFlow/SinglePhaseWell.hpp
     fluidFlow/WellSolverBase.hpp
     multiphysics/HydrofractureSolver.hpp
     multiphysics/PoroelasticSolver.hpp
     multiphysics/ReservoirSolver.hpp
     simplePDE/LaplaceFEM.hpp
     solidMechanics/SolidMechanicsLagrangianFEM.hpp
     solidMechanics/SolidMechanicsLagrangianSSLE.hpp
     solidMechanics/SolidMechanicsLagrangianFEMKernels.hpp
     solidMechanics/SolidMechanicsLagrangianSSLEKernels.hpp
     surfaceGeneration/SurfaceGenerator.hpp
     )

#
# Specify solver sources
#
set( physicsSolvers_sources
     PhysicsSolverManager.cpp
     SolverBase.cpp
     SystemSolverParameters.cpp
     NonlinearSolverParameters.cpp
     fluidFlow/CompositionalMultiphaseFlow.cpp
     fluidFlow/CompositionalMultiphaseFlowKernels.cpp
     fluidFlow/CompositionalMultiphaseWell.cpp
     fluidFlow/FlowSolverBase.cpp
     fluidFlow/SinglePhaseBase.cpp
     fluidFlow/SinglePhaseHybridFVM.cpp
     fluidFlow/SinglePhaseFVM.cpp     
     fluidFlow/SinglePhaseKernels.cpp
     fluidFlow/SinglePhaseWell.cpp
     fluidFlow/WellSolverBase.cpp
     multiphysics/HydrofractureSolver.cpp
     multiphysics/PoroelasticSolver.cpp
     multiphysics/ReservoirSolver.cpp
     simplePDE/LaplaceFEM.cpp
     solidMechanics/SolidMechanicsLagrangianFEM.cpp
     solidMechanics/SolidMechanicsLagrangianSSLE.cpp
     surfaceGeneration/SurfaceGenerator.cpp
     )

if( BUILD_OBJ_LIBS)
  set( dependencyList common constitutive dataRepository linearAlgebra )
else()
  set( dependencyList common linearAlgebra )
endif()


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
                 OBJECT                ${buildAsObj}
               )

target_include_directories( physicsSolvers PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )
if( externalComponentDeps )
  target_include_directories( physicsSolvers PUBLIC ${CMAKE_SOURCE_DIR}/externalComponents)
endif()


geosx_add_code_checks( PREFIX physicsSolvers )

add_subdirectory( fluidFlow/unitTests )

message(STATUS "Leaving src/coreComponents/physicsSolvers/CMakeLists.txt")
