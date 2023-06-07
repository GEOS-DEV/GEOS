#
# Specify solver headers
#
set( physicsSolvers_headers
     LinearSolverParameters.hpp
     NonlinearSolverParameters.hpp
     PhysicsSolverManager.hpp
     SolverBase.hpp
     SolverBaseKernels.hpp
     SolverStatistics.hpp
     FieldStatisticsBase.hpp
     contact/ContactSolverBase.hpp
     contact/ContactFields.hpp
     contact/LagrangianContactSolver.hpp
     contact/SolidMechanicsEFEMKernelsBase.hpp
     contact/SolidMechanicsEFEMKernels.hpp
     contact/SolidMechanicsEFEMStaticCondensationKernels.hpp
     contact/SolidMechanicsEFEMKernelsHelper.hpp
     contact/SolidMechanicsEmbeddedFractures.hpp
     fluidFlow/FlowSolverBase.hpp
     fluidFlow/SinglePhaseBase.hpp
     fluidFlow/SinglePhaseStatistics.hpp
     fluidFlow/SinglePhaseFVM.hpp
     fluidFlow/SinglePhaseHybridFVM.hpp
     fluidFlow/SinglePhaseProppantBase.hpp
     fluidFlow/CompositionalMultiphaseBase.hpp
     fluidFlow/CompositionalMultiphaseFVM.hpp
     fluidFlow/CompositionalMultiphaseHybridFVM.hpp
     fluidFlow/CompositionalMultiphaseUtilities.hpp
     fluidFlow/CompositionalMultiphaseStatistics.hpp
     fluidFlow/ReactiveCompositionalMultiphaseOBL.hpp
     fluidFlow/StencilAccessors.hpp
     fluidFlow/fields/CompositionalMultiphaseBaseFields.hpp
     fluidFlow/fields/ReactiveCompositionalMultiphaseOBLFields.hpp
     fluidFlow/fields/FlowSolverBaseFields.hpp
     fluidFlow/fields/SinglePhaseBaseFields.hpp
     fluidFlow/kernels/MinimumPoreVolumeKernel.hpp
     fluidFlow/kernels/compositionalMultiphase/CapillaryPressureUpdateKernel.hpp
     fluidFlow/kernels/compositionalMultiphase/ComponentFractionKernel.hpp
     fluidFlow/kernels/compositionalMultiphase/ElementBasedAssemblyKernel.hpp
     fluidFlow/kernels/compositionalMultiphase/HydrostaticPressureKernel.hpp
     fluidFlow/kernels/compositionalMultiphase/KernelUtilities.hpp
     fluidFlow/kernels/compositionalMultiphase/PhaseVolumeFractionKernel.hpp
     fluidFlow/kernels/compositionalMultiphase/PropertyKernelBase.hpp
     fluidFlow/kernels/compositionalMultiphase/RelativePermeabilityUpdateKernel.hpp
     fluidFlow/kernels/compositionalMultiphase/ResidualNormKernel.hpp
     fluidFlow/kernels/compositionalMultiphase/ScalingAndCheckingSystemSolutionKernelBase.hpp
     fluidFlow/kernels/compositionalMultiphase/ScalingForSystemSolutionKernel.hpp
     fluidFlow/kernels/compositionalMultiphase/SolutionCheckKernel.hpp
     fluidFlow/kernels/compositionalMultiphase/StatisticsKernel.hpp
     fluidFlow/kernels/compositionalMultiphase/FVM/AquiferBCKernel.hpp
     fluidFlow/kernels/compositionalMultiphase/FVM/CFLKernels.hpp
     fluidFlow/kernels/compositionalMultiphase/FVM/DirichletFaceBasedAssemblyKernel.hpp
     fluidFlow/kernels/compositionalMultiphase/FVM/FaceBasedAssemblyKernel.hpp
     fluidFlow/kernels/compositionalMultiphase/FVM/FaceBasedAssemblyKernelBase.hpp
     fluidFlow/kernels/compositionalMultiphase/FVM/FluxUtilities.hpp
     fluidFlow/kernels/compositionalMultiphase/FVM/PhaseMobilityKernel.hpp
     fluidFlow/kernels/compositionalMultiphase/FVM/stabilized/FaceBasedAssemblyKernel.hpp
     fluidFlow/kernels/compositionalMultiphase/hybridFVM/AssemblerKernel.hpp
     fluidFlow/kernels/compositionalMultiphase/hybridFVM/AssemblerKernelHelper.hpp
     fluidFlow/kernels/compositionalMultiphase/hybridFVM/CellConnectivity.hpp
     fluidFlow/kernels/compositionalMultiphase/hybridFVM/FluxKernel.hpp
     fluidFlow/kernels/compositionalMultiphase/hybridFVM/KernelUtilities.hpp
     fluidFlow/kernels/compositionalMultiphase/hybridFVM/PhaseMobilityKernel.hpp
     fluidFlow/kernels/compositionalMultiphase/hybridFVM/PrecomputeKernel.hpp
     fluidFlow/kernels/compositionalMultiphase/hybridFVM/ResidualNormKernel.hpp
     fluidFlow/kernels/compositionalMultiphase/hybridFVM/SolutionCheckKernel.hpp
     fluidFlow/kernels/compositionalMultiphase/hybridFVM/UpwindingHelper.hpp
     fluidFlow/kernels/compositionalMultiphase/thermal/DirichletFaceBasedAssemblyKernel.hpp
     fluidFlow/kernels/compositionalMultiphase/thermal/ElementBasedAssemblyKernel.hpp
     fluidFlow/kernels/compositionalMultiphase/thermal/FaceBasedAssemblyKernel.hpp
     fluidFlow/kernels/compositionalMultiphase/thermal/FluidUpdateKernel.hpp
     fluidFlow/kernels/compositionalMultiphase/thermal/PhaseMobilityKernel.hpp
     fluidFlow/kernels/compositionalMultiphase/thermal/PhaseVolumeFractionKernel.hpp
     fluidFlow/kernels/compositionalMultiphase/thermal/ResidualNormKernel.hpp
     fluidFlow/kernels/compositionalMultiphase/thermal/ScalingForSystemSolutionKernel.hpp
     fluidFlow/kernels/compositionalMultiphase/thermal/SolidInternalEnergyUpdateKernel.hpp
     fluidFlow/kernels/compositionalMultiphase/thermal/SolutionCheckKernel.hpp
     fluidFlow/kernels/singlePhase/HydrostaticPressureKernel.hpp
     fluidFlow/kernels/singlePhase/SolutionCheckKernel.hpp
     fluidFlow/kernels/singlePhase/StatisticsKernel.hpp
     fluidFlow/kernels/singlePhase/FVM/AquiferBCKernel.hpp
     fluidFlow/kernels/singlePhase/FVM/DirichletFaceBasedAssemblyKernel.hpp
     fluidFlow/kernels/singlePhase/FVM/ElementBasedAssemblyKernel.hpp
     fluidFlow/kernels/singlePhase/FVM/FaceBasedAssemblyKernel.hpp
     fluidFlow/kernels/singlePhase/FVM/FaceBasedAssemblyKernelBase.hpp
     fluidFlow/kernels/singlePhase/FVM/FluxKernel.hpp
     fluidFlow/kernels/singlePhase/FVM/FluxKernelsHelper.hpp
     fluidFlow/kernels/singlePhase/FVM/MobilityKernel.hpp
     fluidFlow/kernels/singlePhase/FVM/ResidualNormKernel.hpp
     fluidFlow/kernels/singlePhase/hybridFVM/ElementBasedAssemblyKernel.hpp
     fluidFlow/kernels/singlePhase/hybridFVM/KernelUtilities.hpp
     fluidFlow/kernels/singlePhase/hybridFVM/ResidualNormKernel.hpp
     fluidFlow/kernels/singlePhase/proppant/FluidUpdateKernel.hpp
     fluidFlow/kernels/singlePhase/thermal/DirichletFaceBasedAssemblyKernel.hpp
     fluidFlow/kernels/singlePhase/thermal/ElementBasedAssemblyKernel.hpp
     fluidFlow/kernels/singlePhase/thermal/FaceBasedAssemblyKernel.hpp
     fluidFlow/kernels/singlePhase/thermal/FluidUpdateKernel.hpp
     fluidFlow/kernels/singlePhase/thermal/MobilityKernel.hpp
     fluidFlow/kernels/singlePhase/thermal/ResidualNormKernel.hpp
     fluidFlow/kernels/singlePhase/thermal/SolidInternalEnergyUpdateKernel.hpp
     fluidFlow/kernels/OBL/ReactiveCompositionalMultiphaseOBLKernels.hpp
     fluidFlow/proppantTransport/ProppantTransport.hpp
     fluidFlow/proppantTransport/ProppantTransportFields.hpp
     fluidFlow/proppantTransport/ProppantTransportKernels.hpp
     fluidFlow/wells/CompositionalMultiphaseWell.hpp
     fluidFlow/wells/CompositionalMultiphaseWellFields.hpp
     fluidFlow/wells/CompositionalMultiphaseWellKernels.hpp
     fluidFlow/wells/SinglePhaseWell.hpp
     fluidFlow/wells/SinglePhaseWellFields.hpp
     fluidFlow/wells/SinglePhaseWellKernels.hpp
     fluidFlow/wells/WellConstants.hpp
     fluidFlow/wells/WellControls.hpp
     fluidFlow/wells/WellSolverBase.hpp
     fluidFlow/wells/WellSolverBaseFields.hpp
     multiphysics/CompositionalMultiphaseReservoirAndWells.hpp
     multiphysics/CoupledReservoirAndWellsBase.hpp
     multiphysics/CoupledSolver.hpp
     multiphysics/FlowProppantTransportSolver.hpp
     multiphysics/HydrofractureSolver.hpp
     multiphysics/HydrofractureSolverKernels.hpp
     multiphysics/MultiphasePoromechanics.hpp
     multiphysics/PhaseFieldFractureSolver.hpp
     multiphysics/PoromechanicsInitialization.hpp
     multiphysics/PoromechanicsFields.hpp
     multiphysics/PoromechanicsInitialization.hpp
     multiphysics/poromechanicsKernels/MultiphasePoromechanics.hpp
     multiphysics/poromechanicsKernels/MultiphasePoromechanics_impl.hpp
     multiphysics/poromechanicsKernels/PoromechanicsBase.hpp
     multiphysics/poromechanicsKernels/SinglePhasePoromechanics.hpp
     multiphysics/poromechanicsKernels/SinglePhasePoromechanics_impl.hpp
     multiphysics/poromechanicsKernels/SinglePhasePoromechanicsEFEM.hpp
     multiphysics/poromechanicsKernels/SinglePhasePoromechanicsEFEM_impl.hpp
     multiphysics/poromechanicsKernels/SinglePhasePoromechanicsFractures.hpp
     multiphysics/poromechanicsKernels/ThermalMultiphasePoromechanics.hpp
     multiphysics/poromechanicsKernels/ThermalMultiphasePoromechanics_impl.hpp
     multiphysics/poromechanicsKernels/ThermalSinglePhasePoromechanics.hpp
     multiphysics/poromechanicsKernels/ThermalSinglePhasePoromechanics_impl.hpp
     multiphysics/SinglePhasePoromechanicsFluxKernels.hpp
     multiphysics/SinglePhasePoromechanics.hpp
     multiphysics/SinglePhasePoromechanicsEmbeddedFractures.hpp
     multiphysics/SinglePhasePoromechanicsConformingFractures.hpp
     multiphysics/SinglePhaseReservoirAndWells.hpp
     simplePDE/LaplaceBaseH1.hpp
     simplePDE/LaplaceFEM.hpp
     simplePDE/LaplaceFEMKernels.hpp
     simplePDE/PhaseFieldDamageFEM.hpp
     simplePDE/PhaseFieldDamageFEMKernels.hpp
     solidMechanics/SolidMechanicsFields.hpp
     solidMechanics/SolidMechanicsLagrangianFEM.hpp
     solidMechanics/SolidMechanicsLagrangianSSLE.hpp
     solidMechanics/kernels/SolidMechanicsLagrangianFEMKernels.hpp
     solidMechanics/kernels/ExplicitFiniteStrain.hpp
     solidMechanics/kernels/ExplicitFiniteStrain_impl.hpp
     solidMechanics/kernels/ExplicitSmallStrain.hpp
     solidMechanics/kernels/ExplicitSmallStrain_impl.hpp
     solidMechanics/kernels/FixedStressThermoPoromechanics.hpp
     solidMechanics/kernels/FixedStressThermoPoromechanics_impl.hpp
     solidMechanics/kernels/ImplicitSmallStrainNewmark.hpp
     solidMechanics/kernels/ImplicitSmallStrainNewmark_impl.hpp
     solidMechanics/kernels/ImplicitSmallStrainQuasiStatic.hpp
     solidMechanics/kernels/ImplicitSmallStrainQuasiStatic_impl.hpp
     solidMechanics/SolidMechanicsStateReset.hpp
     solidMechanics/SolidMechanicsStatistics.hpp
     surfaceGeneration/EmbeddedSurfaceGenerator.hpp
     surfaceGeneration/EmbeddedSurfacesParallelSynchronization.hpp
     surfaceGeneration/ParallelTopologyChange.hpp
     surfaceGeneration/SurfaceGenerator.hpp
     surfaceGeneration/SurfaceGeneratorFields.hpp
     wavePropagation/WaveSolverBase.hpp
     wavePropagation/WaveSolverUtils.hpp
     wavePropagation/WaveSolverBaseFields.hpp
     wavePropagation/AcousticWaveEquationSEM.hpp
     wavePropagation/AcousticWaveEquationSEMKernel.hpp
     wavePropagation/ElasticWaveEquationSEM.hpp
     wavePropagation/ElasticWaveEquationSEMKernel.hpp
     wavePropagation/ElasticFirstOrderWaveEquationSEM.hpp
     wavePropagation/ElasticFirstOrderWaveEquationSEMKernel.hpp
     wavePropagation/AcousticFirstOrderWaveEquationSEM.hpp
     wavePropagation/AcousticFirstOrderWaveEquationSEMKernel.hpp
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
     fluidFlow/FlowSolverBase.cpp
     fluidFlow/SinglePhaseBase.cpp
     fluidFlow/SinglePhaseStatistics.cpp
     fluidFlow/SinglePhaseFVM.cpp
     fluidFlow/SinglePhaseHybridFVM.cpp
     fluidFlow/SinglePhaseProppantBase.cpp
     fluidFlow/CompositionalMultiphaseBase.cpp
     fluidFlow/CompositionalMultiphaseFVM.cpp
     fluidFlow/CompositionalMultiphaseHybridFVM.cpp
     fluidFlow/CompositionalMultiphaseStatistics.cpp
     fluidFlow/ReactiveCompositionalMultiphaseOBL.cpp
     fluidFlow/kernels/compositionalMultiphase/FVM/AquiferBCKernel.cpp
     fluidFlow/kernels/compositionalMultiphase/FVM/CFLKernels.cpp
     fluidFlow/kernels/compositionalMultiphase/FVM/FaceBasedAssemblyKernelBase.cpp
     fluidFlow/kernels/compositionalMultiphase/hybridFVM/AssemblerKernel.cpp
     fluidFlow/kernels/compositionalMultiphase/hybridFVM/AssemblerKernelHelper.cpp
     fluidFlow/kernels/compositionalMultiphase/hybridFVM/FluxKernel.cpp
     fluidFlow/kernels/compositionalMultiphase/hybridFVM/UpwindingHelper.cpp
     fluidFlow/proppantTransport/ProppantTransport.cpp
     fluidFlow/proppantTransport/ProppantTransportKernels.cpp
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
     multiphysics/MultiphasePoromechanics.cpp
     multiphysics/PhaseFieldFractureSolver.cpp
     multiphysics/PoromechanicsInitialization.cpp
     multiphysics/SinglePhasePoromechanics.cpp
     multiphysics/SinglePhasePoromechanicsEmbeddedFractures.cpp
     multiphysics/SinglePhasePoromechanicsConformingFractures.cpp
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
     wavePropagation/ElasticWaveEquationSEM.cpp
     wavePropagation/ElasticFirstOrderWaveEquationSEM.cpp
     wavePropagation/AcousticFirstOrderWaveEquationSEM.cpp

   )

  include( solidMechanics/kernels/SolidMechanicsKernels.cmake)
  include( multiphysics/poromechanicsKernels/PoromechanicsKernels.cmake)

set( dependencyList ${parallelDeps} constitutive mesh linearAlgebra discretizationMethods events )
if( ENABLE_PYGEOSX )
  list( APPEND physicsSolvers_headers
  python/PySolver.hpp
  python/PySolverType.hpp )
  list( APPEND physicsSolvers_sources
        python/PySolver.cpp )
  list( APPEND dependencyList Python3::Python pylvarray )
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

