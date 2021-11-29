#
# Specify solver headers
#
set( physicsSolvers_headers
     LinearSolverParameters.hpp
     NonlinearSolverParameters.hpp
     PhysicsSolverManager.hpp
     SolverBase.hpp
     fluidFlow/CompositionalMultiphaseBase.hpp
     fluidFlow/CompositionalMultiphaseBaseKernels.hpp
     fluidFlow/CompositionalMultiphaseFVM.hpp
     fluidFlow/CompositionalMultiphaseFVMKernels.hpp
     fluidFlow/CompositionalMultiphaseHybridFVM.hpp
     fluidFlow/CompositionalMultiphaseHybridFVMKernels.hpp
     fluidFlow/CompositionalMultiphaseUtilities.hpp
     fluidFlow/FlowSolverBase.hpp
     fluidFlow/FluxKernelsHelper.hpp
     fluidFlow/HybridFVMHelperKernels.hpp          
     fluidFlow/proppantTransport/ProppantTransport.hpp
     fluidFlow/proppantTransport/ProppantTransportKernels.hpp
     fluidFlow/SinglePhaseBase.hpp
     fluidFlow/SinglePhaseBaseKernels.hpp
     fluidFlow/SinglePhaseFVM.hpp
     fluidFlow/SinglePhaseFVMKernels.hpp          
     fluidFlow/SinglePhaseHybridFVM.hpp     
     fluidFlow/SinglePhaseHybridFVMKernels.hpp     
     fluidFlow/SinglePhaseProppantBase.hpp
     fluidFlow/SinglePhaseProppantBaseKernels.hpp
     fluidFlow/wells/CompositionalMultiphaseWell.hpp
     fluidFlow/wells/CompositionalMultiphaseWellKernels.hpp
     fluidFlow/wells/SinglePhaseWell.hpp
     fluidFlow/wells/SinglePhaseWellKernels.hpp     
     fluidFlow/wells/WellControls.hpp
     fluidFlow/wells/WellSolverBase.hpp
     multiphysics/CompositionalMultiphaseReservoir.hpp
     multiphysics/FlowProppantTransportSolver.hpp
     multiphysics/HydrofractureSolver.hpp
     multiphysics/HydrofractureSolverKernels.hpp
     multiphysics/LagrangianContactSolver.hpp
     multiphysics/MultiphasePoromechanicsKernel.hpp
     multiphysics/MultiphasePoromechanicsSolver.hpp
     multiphysics/PhaseFieldFractureSolver.hpp
     multiphysics/ReservoirSolverBase.hpp
     multiphysics/SinglePhasePoromechanicsKernel.hpp
     multiphysics/SinglePhasePoromechanicsEFEMKernel.hpp
     multiphysics/SinglePhasePoromechanicsFluxKernels.hpp
     multiphysics/SinglePhasePoromechanicsSolver.hpp
     multiphysics/SinglePhasePoromechanicsSolverEmbeddedFractures.hpp
     multiphysics/SinglePhaseReservoir.hpp
     simplePDE/LaplaceBaseH1.hpp
     simplePDE/LaplaceFEM.hpp
     simplePDE/LaplaceFEMKernels.hpp
     simplePDE/LaplaceVEM.hpp
     simplePDE/PhaseFieldDamageFEM.hpp
     simplePDE/PhaseFieldDamageFEMKernels.hpp
     solidMechanics/SolidMechanicsEFEMKernels.hpp
     solidMechanics/SolidMechanicsEFEMKernelsHelper.hpp
     solidMechanics/SolidMechanicsEmbeddedFractures.hpp
     solidMechanics/SolidMechanicsFiniteStrainExplicitNewmarkKernel.hpp
     solidMechanics/SolidMechanicsLagrangianFEM.hpp
     solidMechanics/SolidMechanicsLagrangianFEMKernels.hpp
     solidMechanics/SolidMechanicsLagrangianSSLE.hpp
     solidMechanics/SolidMechanicsSmallStrainExplicitNewmarkKernel.hpp
     solidMechanics/SolidMechanicsSmallStrainImplicitNewmarkKernel.hpp
     solidMechanics/SolidMechanicsSmallStrainQuasiStaticKernel.hpp
     surfaceGeneration/EmbeddedSurfaceGenerator.hpp
     surfaceGeneration/EmbeddedSurfacesParallelSynchronization.hpp
     surfaceGeneration/ParallelTopologyChange.hpp
     surfaceGeneration/SurfaceGenerator.hpp
     wavePropagation/WaveSolverBase.hpp
     wavePropagation/AcousticWaveEquationSEM.hpp
     wavePropagation/AcousticWaveEquationSEMKernel.hpp
   )

#
# Specify solver sources
#
set( physicsSolvers_sources
     LinearSolverParameters.cpp
     NonlinearSolverParameters.cpp
     PhysicsSolverManager.cpp
     SolverBase.cpp
     fluidFlow/CompositionalMultiphaseBase.cpp
     fluidFlow/CompositionalMultiphaseBaseKernels.cpp
     fluidFlow/CompositionalMultiphaseFVM.cpp
     fluidFlow/CompositionalMultiphaseFVMKernels.cpp
     fluidFlow/CompositionalMultiphaseHybridFVM.cpp
     fluidFlow/CompositionalMultiphaseHybridFVMKernels.cpp
     fluidFlow/FlowSolverBase.cpp
     fluidFlow/proppantTransport/ProppantTransport.cpp
     fluidFlow/proppantTransport/ProppantTransportKernels.cpp
     fluidFlow/SinglePhaseBase.cpp
     fluidFlow/SinglePhaseBaseKernels.cpp
     fluidFlow/SinglePhaseFVM.cpp
     fluidFlow/SinglePhaseHybridFVM.cpp     
     fluidFlow/SinglePhaseProppantBase.cpp     
     fluidFlow/wells/CompositionalMultiphaseWell.cpp
     fluidFlow/wells/CompositionalMultiphaseWellKernels.cpp     
     fluidFlow/wells/SinglePhaseWell.cpp
     fluidFlow/wells/SinglePhaseWellKernels.cpp
     fluidFlow/wells/WellControls.cpp
     fluidFlow/wells/WellSolverBase.cpp
     multiphysics/CompositionalMultiphaseReservoir.cpp
     multiphysics/FlowProppantTransportSolver.cpp
     multiphysics/HydrofractureSolver.cpp
     multiphysics/LagrangianContactSolver.cpp
     multiphysics/MultiphasePoromechanicsSolver.cpp
     multiphysics/PhaseFieldFractureSolver.cpp
     multiphysics/ReservoirSolverBase.cpp
     multiphysics/SinglePhasePoromechanicsSolver.cpp
     multiphysics/SinglePhasePoromechanicsSolverEmbeddedFractures.cpp
     multiphysics/SinglePhaseReservoir.cpp
     multiphysics/SinglePhasePoromechanicsFluxKernels.cpp
     simplePDE/LaplaceBaseH1.cpp
     simplePDE/LaplaceFEM.cpp
     simplePDE/LaplaceVEM.cpp
     simplePDE/PhaseFieldDamageFEM.cpp
     solidMechanics/SolidMechanicsEmbeddedFractures.cpp
     solidMechanics/SolidMechanicsLagrangianFEM.cpp
     solidMechanics/SolidMechanicsLagrangianSSLE.cpp
     surfaceGeneration/EmbeddedSurfaceGenerator.cpp
     surfaceGeneration/EmbeddedSurfacesParallelSynchronization.cpp
     surfaceGeneration/ParallelTopologyChange.cpp
     surfaceGeneration/SurfaceGenerator.cpp
     wavePropagation/WaveSolverBase.cpp
     wavePropagation/AcousticWaveEquationSEM.cpp
   )

set( dependencyList constitutive mesh linearAlgebra discretizationMethods events )

if ( ENABLE_CUDA )
  set( dependencyList ${dependencyList} cuda )
endif()

blt_add_library( NAME       physicsSolvers
                 SOURCES    ${physicsSolvers_sources}
                 HEADERS    ${physicsSolvers_headers}
                 DEPENDS_ON ${dependencyList} ${externalComponentDeps}
                 OBJECT     ${GEOSX_BUILD_OBJ_LIBS}
               )

target_include_directories( physicsSolvers PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )
if( externalComponentDeps )
  target_include_directories( physicsSolvers PUBLIC ${CMAKE_SOURCE_DIR}/externalComponents )
endif()

geosx_add_code_checks( PREFIX physicsSolvers )
