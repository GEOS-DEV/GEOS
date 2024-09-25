# Specify solver headers
set( physicsSolvers_headers
     ${physicsSolvers_headers}
     fluidFlow/CompositionalMultiphaseBase.hpp
     fluidFlow/CompositionalMultiphaseBaseFields.hpp
     fluidFlow/CompositionalMultiphaseStatistics.hpp
     fluidFlow/kernels/compositional/IsothermalCompositionalMultiphaseBaseKernels.hpp
     fluidFlow/kernels/compositional/ThermalCompositionalMultiphaseBaseKernels.hpp
     fluidFlow/CompositionalMultiphaseFVM.hpp
     fluidFlow/kernels/compositional/IsothermalCompositionalMultiphaseFVMKernels.hpp
     fluidFlow/kernels/compositional/IsothermalCompositionalMultiphaseFVMKernelUtilities.hpp
     fluidFlow/kernels/compositional/ThermalCompositionalMultiphaseFVMKernels.hpp
     fluidFlow/CompositionalMultiphaseHybridFVM.hpp
     fluidFlow/kernels/compositional/CompositionalMultiphaseHybridFVMKernels.hpp
     fluidFlow/CompositionalMultiphaseUtilities.hpp
     fluidFlow/ReactiveCompositionalMultiphaseOBL.hpp
     fluidFlow/ReactiveCompositionalMultiphaseOBLFields.hpp
     fluidFlow/kernels/compositional/ReactiveCompositionalMultiphaseOBLKernels.hpp
     fluidFlow/FlowSolverBase.hpp
     fluidFlow/FlowSolverBaseFields.hpp
     fluidFlow/kernels/FlowSolverBaseKernels.hpp
     fluidFlow/kernels/singlePhase/FluxKernelsHelper.hpp
     fluidFlow/kernels/HybridFVMHelperKernels.hpp
     fluidFlow/SourceFluxStatistics.hpp
     fluidFlow/proppantTransport/ProppantTransport.hpp
     fluidFlow/proppantTransport/ProppantTransportFields.hpp
     fluidFlow/proppantTransport/ProppantTransportKernels.hpp
     fluidFlow/SinglePhaseBase.hpp
     fluidFlow/SinglePhaseBaseFields.hpp
     fluidFlow/kernels/singlePhase/SinglePhaseBaseKernels.hpp
     fluidFlow/SinglePhaseStatistics.hpp
     fluidFlow/SinglePhaseFVM.hpp
     fluidFlow/kernels/singlePhase/SinglePhaseFVMKernels.hpp
     fluidFlow/SinglePhaseHybridFVM.hpp
     fluidFlow/kernels/singlePhase/SinglePhaseHybridFVMKernels.hpp
     fluidFlow/SinglePhaseProppantBase.hpp
     fluidFlow/kernels/singlePhase/SinglePhaseProppantBaseKernels.hpp
     fluidFlow/kernels/singlePhase/SinglePhaseProppantFluxKernels.hpp
     fluidFlow/kernels/compositional/StabilizedCompositionalMultiphaseFVMKernels.hpp
     fluidFlow/kernels/singlePhase/StabilizedSinglePhaseFVMKernels.hpp
     fluidFlow/StencilAccessors.hpp
     fluidFlow/StencilDataCollection.hpp
     fluidFlow/kernels/singlePhase/ThermalSinglePhaseBaseKernels.hpp
     fluidFlow/kernels/singlePhase/ThermalSinglePhaseFVMKernels.hpp
     fluidFlow/wells/CompositionalMultiphaseWell.hpp
     fluidFlow/wells/CompositionalMultiphaseWellFields.hpp
     fluidFlow/wells/kernels/CompositionalMultiphaseWellKernels.hpp
     fluidFlow/wells/SinglePhaseWell.hpp
     fluidFlow/wells/SinglePhaseWellFields.hpp
     fluidFlow/wells/kernels/SinglePhaseWellKernels.hpp
     fluidFlow/wells/WellConstants.hpp
     fluidFlow/wells/WellControls.hpp
     fluidFlow/wells/WellSolverBase.hpp
     fluidFlow/wells/WellSolverBaseFields.hpp
     PARENT_SCOPE )

# Specify solver sources
set( physicsSolvers_sources
     ${physicsSolvers_sources}
     fluidFlow/CompositionalMultiphaseBase.cpp
     fluidFlow/CompositionalMultiphaseFVM.cpp
     fluidFlow/CompositionalMultiphaseStatistics.cpp
     fluidFlow/kernels/compositional/IsothermalCompositionalMultiphaseFVMKernels.cpp
     fluidFlow/CompositionalMultiphaseHybridFVM.cpp
     fluidFlow/kernels/compositional/CompositionalMultiphaseHybridFVMKernels.cpp
     fluidFlow/ReactiveCompositionalMultiphaseOBL.cpp
     fluidFlow/FlowSolverBase.cpp
     fluidFlow/proppantTransport/ProppantTransport.cpp
     fluidFlow/proppantTransport/ProppantTransportKernels.cpp
     fluidFlow/SinglePhaseBase.cpp
     fluidFlow/SinglePhaseStatistics.cpp
     fluidFlow/SinglePhaseFVM.cpp
     fluidFlow/SinglePhaseHybridFVM.cpp
     fluidFlow/SinglePhaseProppantBase.cpp
     fluidFlow/kernels/singlePhase/SinglePhaseProppantFluxKernels.cpp
     fluidFlow/SourceFluxStatistics.cpp
     fluidFlow/StencilDataCollection.cpp
     fluidFlow/wells/CompositionalMultiphaseWell.cpp
     fluidFlow/wells/kernels/CompositionalMultiphaseWellKernels.cpp
     fluidFlow/wells/SinglePhaseWell.cpp
     fluidFlow/wells/kernels/SinglePhaseWellKernels.cpp
     fluidFlow/wells/WellControls.cpp
     fluidFlow/wells/WellSolverBase.cpp 
     PARENT_SCOPE )


