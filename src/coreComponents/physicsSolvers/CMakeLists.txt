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
     multiphysics/MultiphasePoromechanicsKernel.hpp
     multiphysics/MultiphasePoromechanicsSolver.hpp
     multiphysics/PhaseFieldFractureSolver.hpp
     multiphysics/SinglePhasePoromechanicsKernel.hpp
     multiphysics/SinglePhasePoromechanicsEFEMKernel.hpp
     multiphysics/SinglePhasePoromechanicsFluxKernels.hpp
     multiphysics/SinglePhasePoromechanicsSolver.hpp
     multiphysics/SinglePhasePoromechanicsSolverEmbeddedFractures.hpp
     multiphysics/SinglePhaseReservoirAndWells.hpp
     simplePDE/LaplaceBaseH1.hpp
     simplePDE/LaplaceFEM.hpp
     simplePDE/LaplaceFEMKernels.hpp
     simplePDE/PhaseFieldDamageFEM.hpp
     simplePDE/PhaseFieldDamageFEMKernels.hpp
     solidMechanics/SolidMechanicsLagrangianFEM.hpp
     solidMechanics/SolidMechanicsLagrangianSSLE.hpp
     solidMechanics/kernels/SolidMechanicsLagrangianFEMKernels.hpp
     solidMechanics/kernels/ExplicitFiniteStrain.hpp
     solidMechanics/kernels/ExplicitFiniteStrain_impl.hpp
     solidMechanics/kernels/ExplicitSmallStrain.hpp
     solidMechanics/kernels/ExplicitSmallStrain_impl.hpp
     solidMechanics/kernels/ImplicitSmallStrainQuasiStatic.hpp
     solidMechanics/kernels/ImplicitSmallStrainQuasiStatic_impl.hpp
#     solidMechanics/kernels/SolidMechanicsSmallStrainImplicitNewmarkKernel.hpp
     solidMechanics/SolidMechanicsStateReset.hpp
     solidMechanics/SolidMechanicsStatistics.hpp       
     surfaceGeneration/EmbeddedSurfaceGenerator.hpp
     surfaceGeneration/EmbeddedSurfacesParallelSynchronization.hpp
     surfaceGeneration/ParallelTopologyChange.hpp
     surfaceGeneration/SurfaceGenerator.hpp
     wavePropagation/WaveSolverBase.hpp
     wavePropagation/AcousticWaveEquationSEM.hpp
     wavePropagation/AcousticWaveEquationSEMKernel.hpp
     wavePropagation/ElasticWaveEquationSEM.hpp
     wavePropagation/ElasticWaveEquationSEMKernel.hpp
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
     multiphysics/MultiphasePoromechanicsSolver.cpp
     multiphysics/PhaseFieldFractureSolver.cpp
     multiphysics/SinglePhasePoromechanicsSolver.cpp
     multiphysics/SinglePhasePoromechanicsSolverEmbeddedFractures.cpp
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
   )

  include( solidMechanics/kernels/SolidMechanicsKernels.cmake)

  # set( kernelNames ExplicitSmallStrain ExplicitFiniteStrain ImplicitSmallStrainQuasiStatic )
  # set( subregionList CellElementSubRegion )
  # set( solidBaseDispatch "DamageSpectral<ElasticIsotropic>;DamageVolDev<ElasticIsotropic>;Damage<ElasticIsotropic>;DruckerPragerExtended;ModifiedCamClay;DelftEgg;DruckerPrager;ElasticIsotropic;ElasticTransverseIsotropic;ElasticIsotropicPressureDependent;ElasticOrthotropic")
  # set( finiteElementDispatch H1_Hexahedron_Lagrange1_GaussLegendre2
  #                            H1_Wedge_Lagrange1_Gauss6
  #                            H1_Tetrahedron_Lagrange1_Gauss1
  #                            H1_Pyramid_Lagrange1_Gauss5
  #                            Q3_Hexahedron_Lagrange_GaussLobatto
  #                            H1_Tetrahedron_VEM_Gauss1
  #                            H1_Wedge_VEM_Gauss1
  #                            H1_Hexahedron_VEM_Gauss1
  #                            H1_Prism5_VEM_Gauss1
  #                            H1_Prism6_VEM_Gauss1
  #                            H1_Prism7_VEM_Gauss1
  #                            H1_Prism8_VEM_Gauss1
  #                            H1_Prism9_VEM_Gauss1
  #                            H1_Prism10_VEM_Gauss1
  #                            H1_Prism11_VEM_Gauss1 )


  # foreach( kernelName ${kernelNames} )
  #   foreach( SUBREGION_TYPE  ${subregionList} )
  #     foreach( CONSTITUTIVE_TYPE ${solidBaseDispatch} )
  #       foreach( FE_TYPE ${finiteElementDispatch} )
  #         message( "SUBREGION_TYPE/CONSTITUTIVE_TYPE/FE_TYPE = ${SUBREGION_TYPE}/${CONSTITUTIVE_TYPE}/${FE_TYPE}")
  #         file(READ "${CMAKE_CURRENT_SOURCE_DIR}/solidMechanics/kernels/SolidMechanicsKernels.cpp.template" FILE_CONTENTS)
  #         string(REPLACE "KERNELNAME" "${kernelName}" FILE_CONTENTS "${FILE_CONTENTS}")
  #         string(REPLACE "SUBREGION_TYPE" "${SUBREGION_TYPE}" FILE_CONTENTS "${FILE_CONTENTS}")
  #         string(REPLACE "CONSTITUTIVE_TYPE" "${CONSTITUTIVE_TYPE}" FILE_CONTENTS "${FILE_CONTENTS}")
  #         string(REPLACE "FE_TYPE" "${FE_TYPE}" FILE_CONTENTS "${FILE_CONTENTS}")
  #         set( filename "solidMechanics/kernels/generatedSrc/${kernelName}_${SUBREGION_TYPE}_${CONSTITUTIVE_TYPE}_${FE_TYPE}.cpp" )
  #         string(REPLACE "<" "" filename "${filename}")
  #         string(REPLACE ">" "" filename "${filename}")
  #         set( filepath "${CMAKE_CURRENT_SOURCE_DIR}/${filename}")
  #         file(WRITE ${filepath} "${FILE_CONTENTS}")
  #         list( APPEND physicsSolvers_sources ${filename} )
  #       endforeach()
  #     endforeach()
  #   endforeach()
  # endforeach()




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

if( ENABLE_CUDA_NVTOOLSEXT )
  set( dependencyList ${dependencyList} CUDA::nvToolsExt )
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

