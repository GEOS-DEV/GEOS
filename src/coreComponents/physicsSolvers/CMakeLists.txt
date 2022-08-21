#
# Specify solver headers
#
set( physicsSolvers_headers
     LinearSolverParameters.hpp
     NonlinearSolverParameters.hpp
     PhysicsSolverManager.hpp
     SolverBase.hpp
     SolverStatistics.hpp     
     FieldStatisticsBase.hpp
     contact/ContactSolverBase.hpp
     contact/ContactExtrinsicData.hpp
     contact/LagrangianContactSolver.hpp
     contact/SolidMechanicsEFEMKernelsBase.hpp
     contact/SolidMechanicsEFEMKernels.hpp
     contact/SolidMechanicsEFEMStaticCondensationKernels.hpp
     contact/SolidMechanicsEFEMKernelsHelper.hpp
     contact/SolidMechanicsEmbeddedFractures.hpp
     fluidFlow/CompositionalMultiphaseBase.hpp
     fluidFlow/CompositionalMultiphaseBaseExtrinsicData.hpp
     fluidFlow/CompositionalMultiphaseStatistics.hpp     
     fluidFlow/IsothermalCompositionalMultiphaseBaseKernels.hpp
     fluidFlow/ThermalCompositionalMultiphaseBaseKernels.hpp     
     fluidFlow/CompositionalMultiphaseFVM.hpp
     fluidFlow/IsothermalCompositionalMultiphaseFVMKernels.hpp
     fluidFlow/ThermalCompositionalMultiphaseFVMKernels.hpp     
     fluidFlow/CompositionalMultiphaseHybridFVM.hpp
     fluidFlow/CompositionalMultiphaseHybridFVMKernels.hpp
     fluidFlow/CompositionalMultiphaseUtilities.hpp
     fluidFlow/ReactiveCompositionalMultiphaseOBL.hpp
     fluidFlow/ReactiveCompositionalMultiphaseOBLExtrinsicData.hpp
     fluidFlow/ReactiveCompositionalMultiphaseOBLKernels.hpp
     fluidFlow/FlowSolverBase.hpp
     fluidFlow/FlowSolverBaseExtrinsicData.hpp
     fluidFlow/FluxKernelsHelper.hpp
     fluidFlow/HybridFVMHelperKernels.hpp
     fluidFlow/proppantTransport/ProppantTransport.hpp
     fluidFlow/proppantTransport/ProppantTransportExtrinsicData.hpp
     fluidFlow/proppantTransport/ProppantTransportKernels.hpp
     fluidFlow/SinglePhaseBase.hpp
     fluidFlow/SinglePhaseBaseExtrinsicData.hpp
     fluidFlow/SinglePhaseBaseKernels.hpp
     fluidFlow/SinglePhaseStatistics.hpp
     fluidFlow/SinglePhaseFVM.hpp
     fluidFlow/SinglePhaseFVMKernels.hpp
     fluidFlow/SinglePhaseHybridFVM.hpp
     fluidFlow/SinglePhaseHybridFVMKernels.hpp
     fluidFlow/SinglePhaseProppantBase.hpp
     fluidFlow/SinglePhaseProppantBaseKernels.hpp
     fluidFlow/StencilAccessors.hpp
     fluidFlow/ThermalSinglePhaseBaseKernels.hpp
     fluidFlow/ThermalSinglePhaseFVMKernels.hpp
     fluidFlow/wells/CompositionalMultiphaseWell.hpp
     fluidFlow/wells/CompositionalMultiphaseWellExtrinsicData.hpp
     fluidFlow/wells/CompositionalMultiphaseWellKernels.hpp
     fluidFlow/wells/SinglePhaseWell.hpp
     fluidFlow/wells/SinglePhaseWellExtrinsicData.hpp
     fluidFlow/wells/SinglePhaseWellKernels.hpp
     fluidFlow/wells/WellControls.hpp
     fluidFlow/wells/WellSolverBase.hpp
     fluidFlow/wells/WellSolverBaseExtrinsicData.hpp
     multiphysics/CompositionalMultiphaseReservoirAndWells.hpp          
     multiphysics/CoupledReservoirAndWellsBase.hpp     
     multiphysics/CoupledSolver.hpp
     multiphysics/FlowProppantTransportSolver.hpp
     multiphysics/HydrofractureSolver.hpp
     multiphysics/HydrofractureSolverKernels.hpp
     multiphysics/LagrangianContactFlowSolver.hpp
     multiphysics/MultiphasePoromechanicsKernel.hpp
     multiphysics/MultiphasePoromechanicsSolver.hpp
     multiphysics/PhaseFieldFractureSolver.hpp
     multiphysics/SinglePhasePoromechanicsKernel.hpp
     multiphysics/SinglePhasePoromechanicsEFEMKernel.hpp
     multiphysics/SinglePhasePoromechanicsFluxKernels.hpp
     multiphysics/SinglePhasePoromechanicsSolver.hpp
     multiphysics/SinglePhasePoromechanicsSolverEmbeddedFractures.hpp
     multiphysics/SinglePhasePoromechanicsLagrangianContactSolver.hpp
     multiphysics/SinglePhaseReservoirAndWells.hpp
     simplePDE/LaplaceBaseH1.hpp
     simplePDE/LaplaceFEM.hpp
     simplePDE/LaplaceFEMKernels.hpp
     simplePDE/PhaseFieldDamageFEM.hpp
     simplePDE/PhaseFieldDamageFEMKernels.hpp
     solidMechanics/SolidMechanicsFiniteStrainExplicitNewmarkKernel.hpp
     solidMechanics/SolidMechanicsLagrangianFEM.hpp
     solidMechanics/SolidMechanicsLagrangianFEMKernels.hpp
     solidMechanics/SolidMechanicsLagrangianSSLE.hpp
     solidMechanics/SolidMechanicsSmallStrainExplicitNewmarkKernel.hpp
     solidMechanics/SolidMechanicsSmallStrainImplicitNewmarkKernel.hpp
     solidMechanics/SolidMechanicsSmallStrainQuasiStaticKernel.hpp
     solidMechanics/SolidMechanicsStateReset.hpp
     solidMechanics/SolidMechanicsStatistics.hpp	     
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
     SolverStatistics.cpp     
     contact/ContactSolverBase.cpp
     contact/LagrangianContactSolver.cpp
     contact/SolidMechanicsEmbeddedFractures.cpp
     fluidFlow/CompositionalMultiphaseBase.cpp
     fluidFlow/CompositionalMultiphaseFVM.cpp
     fluidFlow/CompositionalMultiphaseStatistics.cpp          
     fluidFlow/IsothermalCompositionalMultiphaseFVMKernels.cpp
     fluidFlow/CompositionalMultiphaseHybridFVM.cpp
     fluidFlow/CompositionalMultiphaseHybridFVMKernels.cpp
     fluidFlow/ReactiveCompositionalMultiphaseOBL.cpp
     fluidFlow/FlowSolverBase.cpp
     fluidFlow/proppantTransport/ProppantTransport.cpp
     fluidFlow/proppantTransport/ProppantTransportKernels.cpp
     fluidFlow/SinglePhaseBase.cpp
     fluidFlow/SinglePhaseStatistics.cpp
     fluidFlow/SinglePhaseFVM.cpp
     fluidFlow/SinglePhaseHybridFVM.cpp
     fluidFlow/SinglePhaseProppantBase.cpp
     fluidFlow/wells/CompositionalMultiphaseWell.cpp
     fluidFlow/wells/CompositionalMultiphaseWellKernels.cpp
     fluidFlow/wells/SinglePhaseWell.cpp
     fluidFlow/wells/SinglePhaseWellKernels.cpp
     fluidFlow/wells/WellControls.cpp
     fluidFlow/wells/WellSolverBase.cpp
     multiphysics/CompositionalMultiphaseReservoirAndWells.cpp               
     multiphysics/CoupledReservoirAndWellsBase.cpp		
     multiphysics/FlowProppantTransportSolver.cpp
     multiphysics/HydrofractureSolver.cpp
     multiphysics/LagrangianContactFlowSolver.cpp
     multiphysics/MultiphasePoromechanicsSolver.cpp
     multiphysics/PhaseFieldFractureSolver.cpp
     multiphysics/SinglePhasePoromechanicsSolver.cpp
     multiphysics/SinglePhasePoromechanicsSolverEmbeddedFractures.cpp
     multiphysics/SinglePhasePoromechanicsLagrangianContactSolver.cpp
     multiphysics/SinglePhasePoromechanicsFluxKernels.cpp
     multiphysics/SinglePhaseReservoirAndWells.cpp
     simplePDE/LaplaceBaseH1.cpp
     simplePDE/LaplaceFEM.cpp
     simplePDE/PhaseFieldDamageFEM.cpp
     solidMechanics/SolidMechanicsLagrangianFEM.cpp
     solidMechanics/SolidMechanicsLagrangianSSLE.cpp
     solidMechanics/SolidMechanicsStateReset.cpp
     solidMechanics/SolidMechanicsStatistics.cpp	          
     surfaceGeneration/EmbeddedSurfaceGenerator.cpp
     surfaceGeneration/EmbeddedSurfacesParallelSynchronization.cpp
     surfaceGeneration/ParallelTopologyChange.cpp
     surfaceGeneration/SurfaceGenerator.cpp
     wavePropagation/WaveSolverBase.cpp
     wavePropagation/AcousticWaveEquationSEM.cpp
   )

set( dependencyList constitutive mesh linearAlgebra discretizationMethods events )
if( ENABLE_PYGEOSX )
  list( APPEND physicsSolvers_headers
	python/PySolver.hpp 
	python/PySolverType.hpp ) 
  list( APPEND physicsSolvers_sources
	python/PySolver.cpp ) 
  list( APPEND dependencyList Python3::Python pylvarray )
endif()

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

