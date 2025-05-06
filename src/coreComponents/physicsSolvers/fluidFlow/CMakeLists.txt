# SPDX-License-Identifier: LGPL-2.1-only
#
# Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
# Copyright (c) 2018-2024 TotalEnergies
# Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
# Copyright (c) 2023-2024 Chevron
# Copyright (c) 2019-     GEOS/GEOSX Contributors
# All rights reserved
#
# See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
#
#--------------------------------------------------------------------------------------------------

#[[
Package: fluidFlowSolvers
Contains:
#]]

# Specify solver headers
set( fluidFlowSolvers_headers
     FlowSolverBase.hpp
     FlowSolverBaseFields.hpp
     CompositionalMultiphaseBase.hpp
     CompositionalMultiphaseBaseFields.hpp
     CompositionalMultiphaseStatistics.hpp
     CompositionalMultiphaseFVM.hpp
     CompositionalMultiphaseHybridFVM.hpp
     CompositionalMultiphaseUtilities.hpp
     ImmiscibleMultiphaseFlow.hpp
     ImmiscibleMultiphaseFlowFields.hpp
     LogLevelsInfo.hpp
     ReactiveCompositionalMultiphaseOBL.hpp
     ReactiveCompositionalMultiphaseOBLFields.hpp
     SourceFluxStatistics.hpp
     SinglePhaseBase.hpp
     SinglePhaseBaseFields.hpp
     SinglePhaseStatistics.hpp
     SinglePhaseFVM.hpp
     SinglePhaseHybridFVM.hpp
     SinglePhaseProppantBase.hpp
     StencilAccessors.hpp
     StencilDataCollection.hpp
     LogLevelsInfo.hpp
     kernels/MinPoreVolumeMaxPorosityKernel.hpp
     kernels/StencilWeightsUpdateKernel.hpp
     kernels/HybridFVMHelperKernels.hpp
     kernels/singlePhase/AccumulationKernels.hpp
     kernels/singlePhase/AquiferBCKernel.hpp
     kernels/singlePhase/DirichletFluxComputeKernel.hpp
     kernels/singlePhase/FluidUpdateKernel.hpp
     kernels/singlePhase/FluxComputeKernel.hpp
     kernels/singlePhase/FluxComputeKernelBase.hpp
     kernels/singlePhase/FluxKernelsHelper.hpp
     kernels/singlePhase/HydrostaticPressureKernel.hpp
     kernels/singlePhase/MobilityKernel.hpp
     kernels/singlePhase/ResidualNormKernel.hpp
     kernels/singlePhase/SinglePhaseHybridFVMKernels.hpp
     kernels/singlePhase/SolidInternalEnergyUpdateKernel.hpp
     kernels/singlePhase/SolutionCheckKernel.hpp
     kernels/singlePhase/SolutionScalingKernel.hpp
     kernels/singlePhase/StabilizedFluxComputeKernel.hpp
     kernels/singlePhase/StatisticsKernel.hpp
     kernels/singlePhase/ThermalAccumulationKernels.hpp
     kernels/singlePhase/ThermalDirichletFluxComputeKernel.hpp
     kernels/singlePhase/ThermalFluxComputeKernel.hpp
     kernels/singlePhase/proppant/ProppantBaseKernels.hpp
     kernels/singlePhase/proppant/ProppantFluxKernels.hpp
     kernels/compositional/AccumulationKernel.hpp
     kernels/compositional/AquiferBCKernel.hpp
     kernels/compositional/PPUPhaseFlux.hpp
     kernels/compositional/C1PPUPhaseFlux.hpp
     kernels/compositional/CapillaryPressureUpdateKernel.hpp
     kernels/compositional/CFLKernel.hpp
     kernels/compositional/CompositionalMultiphaseHybridFVMKernels.hpp
     kernels/compositional/DiffusionDispersionFluxComputeKernel.hpp
     kernels/compositional/DirichletFluxComputeKernel.hpp
     kernels/compositional/DissipationFluxComputeKernel.hpp
     kernels/compositional/FluidUpdateKernel.hpp
     kernels/compositional/FluidUpdateKernel.hpp
     kernels/compositional/FluxComputeKernel.hpp
     kernels/compositional/FluxComputeKernelBase.hpp
     kernels/compositional/GlobalComponentFractionKernel.hpp
     kernels/compositional/HydrostaticPressureKernel.hpp
     kernels/compositional/IHUPhaseFlux.hpp
     kernels/compositional/HU2PhaseFlux.hpp
     kernels/compositional/KernelLaunchSelectors.hpp
     kernels/compositional/PhaseComponentFlux.hpp
     kernels/compositional/PhaseMobilityKernel.hpp
     kernels/compositional/PhaseVolumeFractionKernel.hpp
     kernels/compositional/PotGrad.hpp
     kernels/compositional/PPUPhaseFlux.hpp
     kernels/compositional/PropertyKernelBase.hpp
     kernels/compositional/ReactiveCompositionalMultiphaseOBLKernels.hpp
     kernels/compositional/RelativePermeabilityUpdateKernel.hpp
     kernels/compositional/ResidualNormKernel.hpp
     kernels/compositional/SolidInternalEnergyUpdateKernel.hpp
     kernels/compositional/SolutionScalingAndCheckingKernelBase.hpp
     kernels/compositional/SolutionCheckKernel.hpp
     kernels/compositional/SolutionScalingKernel.hpp
     kernels/compositional/StabilizedFluxComputeKernel.hpp
     kernels/compositional/StatisticsKernel.hpp
     kernels/compositional/ThermalAccumulationKernel.hpp
     kernels/compositional/ThermalDiffusionDispersionFluxComputeKernel.hpp
     kernels/compositional/ThermalDirichletFluxComputeKernel.hpp
     kernels/compositional/ThermalFluxComputeKernel.hpp
     kernels/compositional/ThermalPhaseMobilityKernel.hpp
     kernels/compositional/ThermalPhaseVolumeFractionKernel.hpp
     kernels/compositional/ThermalResidualNormKernel.hpp
     kernels/compositional/ThermalSolutionCheckKernel.hpp
     kernels/compositional/ThermalSolutionScalingKernel.hpp
     kernels/immiscibleMultiphase/ImmiscibleMultiphaseKernels.hpp
     wells/CompositionalMultiphaseWell.hpp
     wells/CompositionalMultiphaseWellFields.hpp
     wells/SinglePhaseWell.hpp
     wells/SinglePhaseWellFields.hpp
     wells/WellConstants.hpp
     wells/WellControls.hpp
     wells/WellSolverBase.hpp
     wells/WellSolverBaseFields.hpp
     wells/LogLevelsInfo.hpp
     wells/kernels/SinglePhaseWellKernels.hpp
     wells/kernels/CompositionalMultiphaseWellKernels.hpp
     proppantTransport/ProppantTransport.hpp
     proppantTransport/ProppantTransportFields.hpp
     proppantTransport/ProppantTransportKernels.hpp )

# Specify solver sources
set( fluidFlowSolvers_sources
     CompositionalMultiphaseBase.cpp
     CompositionalMultiphaseFVM.cpp
     CompositionalMultiphaseStatistics.cpp
     CompositionalMultiphaseHybridFVM.cpp
     ImmiscibleMultiphaseFlow.cpp
     ReactiveCompositionalMultiphaseOBL.cpp
     FlowSolverBase.cpp
     SinglePhaseBase.cpp
     SinglePhaseStatistics.cpp
     SinglePhaseFVM.cpp
     SinglePhaseHybridFVM.cpp
     SinglePhaseProppantBase.cpp
     SourceFluxStatistics.cpp
     StencilDataCollection.cpp
     kernels/singlePhase/proppant/ProppantFluxKernels.cpp
     kernels/compositional/AquiferBCKernel.cpp
     kernels/compositional/CFLKernel.cpp
     kernels/compositional/FluxComputeKernelBase.cpp
     kernels/compositional/CompositionalMultiphaseHybridFVMKernels.cpp
     wells/CompositionalMultiphaseWell.cpp
     wells/kernels/CompositionalMultiphaseWellKernels.cpp
     wells/SinglePhaseWell.cpp
     wells/kernels/SinglePhaseWellKernels.cpp
     wells/WellControls.cpp
     wells/WellSolverBase.cpp 
     proppantTransport/ProppantTransport.cpp
     proppantTransport/ProppantTransportKernels.cpp )

set( dependencyList ${parallelDeps} physicsSolversBase )

geos_decorate_link_dependencies( LIST decoratedDependencies DEPENDENCIES ${dependencyList} )




file( READ "${CMAKE_CURRENT_LIST_DIR}/kernelSpecs.json" kernelSpecs )
set( kernelTemplateFileList "" )

list( APPEND kernelTemplateFileList
      SinglePhaseHybridFVMKernels.cpp.template
      CompositionalMultiphaseHybridFVMKernels_upwinding.cpp.template
      CompositionalMultiphaseHybridFVMKernels_assembly.cpp.template
      CompositionalMultiphaseHybridFVMKernels_flux.cpp.template
      ReactiveCompositionalMultiphaseOBLKernels.cpp.template )


foreach( kernelTemplateFile ${kernelTemplateFileList} )
  get_filename_component( jsonKey ${kernelTemplateFile} NAME_WE )
  generateKernels( TEMPLATE ${kernelTemplateFile}
                   JSON kernelSpecs
                   KEY ${jsonKey}
                   HEADERS headerFiles
                   SOURCES sourceFiles )

list(APPEND fluidFlowSolvers_headers ${headerFiles})
list(APPEND fluidFlowSolvers_sources ${sourceFiles})
endforeach()

# TODO: The two kernels below have non-matching file names and JSON keys.
#       Either fix them to follow pattern, or come up with another mechanism.

generateKernels( TEMPLATE "CompositionalMultiphaseFVMKernel_isothermal.cpp.template"
                 JSON kernelSpecs
                 KEY "IsothermalCompositionalMultiphaseFVM_DirichletFluxComputeKernels"
                 HEADERS headerFiles
                 SOURCES sourceFiles )
list(APPEND fluidFlowSolvers_headers ${headerFiles})
list(APPEND fluidFlowSolvers_sources ${sourceFiles})

generateKernels( TEMPLATE "CompositionalMultiphaseFVMKernel_thermal.cpp.template"
                 JSON kernelSpecs
                 KEY "ThermalCompositionalMultiphaseFVM_DirichletFluxComputeKernels"
                 HEADERS headerFiles
                 SOURCES sourceFiles )
list(APPEND fluidFlowSolvers_headers ${headerFiles})
list(APPEND fluidFlowSolvers_sources ${sourceFiles})





blt_add_library( NAME       fluidFlowSolvers
                 SOURCES    ${fluidFlowSolvers_sources}
                 HEADERS    ${fluidFlowSolvers_headers}
                 DEPENDS_ON ${decoratedDependencies} ${externalComponentDeps}
                 OBJECT     ${GEOS_BUILD_OBJ_LIBS}
                 SHARED     ${GEOS_BUILD_SHARED_LIBS}
               )

target_include_directories( fluidFlowSolvers PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )

install( TARGETS fluidFlowSolvers LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib )

if( externalComponentDeps )
  target_include_directories( fluidFlowSolvers PUBLIC ${CMAKE_SOURCE_DIR}/externalComponents )
endif()
