#
# Specify all headers
#
set( constitutive_headers
     ConstitutiveBase.hpp
     ConstitutiveManager.hpp
     ConstitutivePassThru.hpp
     ConstitutivePassThruHandler.hpp
     ExponentialRelation.hpp
     NullModel.hpp
     capillaryPressure/BrooksCoreyCapillaryPressure.hpp
     capillaryPressure/CapillaryPressureBase.hpp
     capillaryPressure/CapillaryPressureExtrinsicData.hpp
     capillaryPressure/JFunctionCapillaryPressure.hpp          
     capillaryPressure/TableCapillaryPressure.hpp
     capillaryPressure/TableCapillaryPressureHelpers.hpp          
     capillaryPressure/VanGenuchtenCapillaryPressure.hpp
     capillaryPressure/capillaryPressureSelector.hpp
     capillaryPressure/layouts.hpp
     contact/ContactBase.hpp
     contact/CoulombContact.hpp
     contact/ContactSelector.hpp
     contact/FrictionlessContact.hpp
     fluid/BlackOilFluidBase.hpp	
     fluid/BlackOilFluid.hpp
     fluid/CompressibleSinglePhaseFluid.hpp
     fluid/DeadOilFluid.hpp
     fluid/MultiFluidBase.hpp
     fluid/MultiFluidUtils.hpp
     fluid/MultiFluidExtrinsicData.hpp
     fluid/MultiPhaseMultiComponentFluid.hpp     
     fluid/PVTDriver.hpp
     fluid/PVTOData.hpp
     fluid/PVTFunctions/PhillipsBrineDensity.hpp
     fluid/PVTFunctions/PhillipsBrineViscosity.hpp
     fluid/PVTFunctions/EzrokhiBrineDensity.hpp
     fluid/PVTFunctions/EzrokhiBrineViscosity.hpp
     fluid/PVTFunctions/CO2Solubility.hpp
     fluid/PVTFunctions/FenghourCO2Viscosity.hpp
     fluid/PVTFunctions/FlashModelBase.hpp
     fluid/PVTFunctions/PVTFunctionBase.hpp
     fluid/PVTFunctions/PVTFunctionHelpers.hpp
     fluid/PVTFunctions/SpanWagnerCO2Density.hpp
     fluid/PVTFunctions/BrineEnthalpy.hpp
     fluid/PVTFunctions/CO2Enthalpy.hpp
     fluid/PVTFunctions/BrineInternalEnergy.hpp
     fluid/PVTFunctions/CO2InternalEnergy.hpp
     fluid/PVTFunctions/CO2EOSSolver.hpp     
     fluid/ParticleFluid.hpp
     fluid/ParticleFluidBase.hpp
     fluid/ProppantSlurryFluid.hpp
     fluid/SingleFluidBase.hpp
     fluid/SingleFluidExtrinsicData.hpp
     fluid/SlurryFluidBase.hpp
     fluid/SlurryFluidExtrinsicData.hpp     
     fluid/layouts.hpp
     fluid/multiFluidSelector.hpp
     fluid/particleFluidSelector.hpp
     fluid/singleFluidSelector.hpp
     fluid/slurryFluidSelector.hpp
     permeability/CarmanKozenyPermeability.hpp
     permeability/ConstantPermeability.hpp
     permeability/ParallelPlatesPermeability.hpp
     permeability/PermeabilityBase.hpp
     permeability/PermeabilityExtrinsicData.hpp     
     permeability/ProppantPermeability.hpp     
     permeability/SlipDependentPermeability.hpp 
     relativePermeability/BrooksCoreyBakerRelativePermeability.hpp
     relativePermeability/BrooksCoreyRelativePermeability.hpp
     relativePermeability/RelativePermeabilityBase.hpp 
     relativePermeability/RelativePermeabilityExtrinsicData.hpp    
     relativePermeability/RelativePermeabilityInterpolators.hpp
     relativePermeability/TableRelativePermeability.hpp
     relativePermeability/VanGenuchtenBakerRelativePermeability.hpp
     relativePermeability/layouts.hpp
     relativePermeability/relativePermeabilitySelector.hpp
     solid/CompressibleSolid.hpp
     solid/ProppantSolid.hpp
     solid/CoupledSolidBase.hpp
     solid/CoupledSolid.hpp
     solid/Damage.hpp
     solid/DamageSpectral.hpp
     solid/DamageSpectralUtilities.hpp
     solid/DamageVolDev.hpp
     solid/DruckerPrager.hpp
     solid/DruckerPragerExtended.hpp
     solid/ModifiedCamClay.hpp
     solid/DelftEgg.hpp
     solid/ElasticIsotropic.hpp
     solid/ElasticIsotropicPressureDependent.hpp
     solid/ElasticTransverseIsotropic.hpp
     solid/ElasticOrthotropic.hpp
     solid/InvariantDecompositions.hpp
     solid/PorousSolid.hpp
     solid/PropertyConversions.hpp
     solid/SolidBase.hpp
     solid/SolidInternalEnergy.hpp
     solid/SolidModelDiscretizationOps.hpp
     solid/SolidModelDiscretizationOpsFullyAnisotroipic.hpp
     solid/SolidModelDiscretizationOpsIsotropic.hpp
     solid/SolidModelDiscretizationOpsTransverseIsotropic.hpp
     solid/SolidModelDiscretizationOpsOrthotropic.hpp
     solid/TriaxialDriver.hpp
     solid/porosity/BiotPorosity.hpp
     solid/porosity/PorosityBase.hpp
     solid/porosity/PressurePorosity.hpp
     solid/porosity/ProppantPorosity.hpp
   )
#
# Specify all sources
#
set( constitutive_sources
     ConstitutiveBase.cpp
     ConstitutiveManager.cpp
     NullModel.cpp
     capillaryPressure/BrooksCoreyCapillaryPressure.cpp
     capillaryPressure/CapillaryPressureBase.cpp
     capillaryPressure/JFunctionCapillaryPressure.cpp               
     capillaryPressure/TableCapillaryPressure.cpp
     capillaryPressure/TableCapillaryPressureHelpers.cpp     
     capillaryPressure/VanGenuchtenCapillaryPressure.cpp
     contact/ContactBase.cpp
     contact/CoulombContact.cpp
     contact/FrictionlessContact.cpp
     fluid/CompressibleSinglePhaseFluid.cpp
     fluid/BlackOilFluidBase.cpp
     fluid/BlackOilFluid.cpp
     fluid/DeadOilFluid.cpp
     fluid/MultiFluidBase.cpp
     fluid/MultiPhaseMultiComponentFluid.cpp
     fluid/PVTDriver.cpp
     fluid/PVTOData.cpp
     fluid/PVTFunctions/PhillipsBrineDensity.cpp
     fluid/PVTFunctions/PhillipsBrineViscosity.cpp
     fluid/PVTFunctions/EzrokhiBrineDensity.cpp
     fluid/PVTFunctions/EzrokhiBrineViscosity.cpp
     fluid/PVTFunctions/CO2Solubility.cpp
     fluid/PVTFunctions/FenghourCO2Viscosity.cpp
     fluid/PVTFunctions/SpanWagnerCO2Density.cpp
     fluid/PVTFunctions/PVTFunctionHelpers.cpp
     fluid/PVTFunctions/BrineEnthalpy.cpp
     fluid/PVTFunctions/CO2Enthalpy.cpp
     fluid/PVTFunctions/BrineInternalEnergy.cpp
     fluid/PVTFunctions/CO2InternalEnergy.cpp
     fluid/PVTFunctions/CO2EOSSolver.cpp          
     fluid/ParticleFluid.cpp              
     fluid/ParticleFluidBase.cpp
     fluid/ProppantSlurryFluid.cpp
     fluid/SingleFluidBase.cpp
     fluid/SlurryFluidBase.cpp
     permeability/CarmanKozenyPermeability.cpp
     permeability/ConstantPermeability.cpp
     permeability/ParallelPlatesPermeability.cpp
     permeability/PermeabilityBase.cpp
     permeability/ProppantPermeability.cpp     
     permeability/SlipDependentPermeability.cpp
     relativePermeability/BrooksCoreyBakerRelativePermeability.cpp
     relativePermeability/BrooksCoreyRelativePermeability.cpp
     relativePermeability/RelativePermeabilityBase.cpp
     relativePermeability/TableRelativePermeability.cpp        
     relativePermeability/VanGenuchtenBakerRelativePermeability.cpp
     solid/CompressibleSolid.cpp
     solid/CoupledSolidBase.cpp
     solid/ProppantSolid.cpp
     solid/Damage.cpp
     solid/DamageSpectral.cpp
     solid/DamageVolDev.cpp
     solid/DruckerPrager.cpp
     solid/DruckerPragerExtended.cpp
     solid/ModifiedCamClay.cpp
     solid/DelftEgg.cpp
     solid/ElasticIsotropic.cpp
     solid/ElasticIsotropicPressureDependent.cpp
     solid/ElasticTransverseIsotropic.cpp
     solid/ElasticOrthotropic.cpp
     solid/PorousSolid.cpp
     solid/SolidBase.cpp
     solid/SolidInternalEnergy.cpp
     solid/TriaxialDriver.cpp
     solid/porosity/BiotPorosity.cpp
     solid/porosity/PorosityBase.cpp
     solid/porosity/PressurePorosity.cpp
     solid/porosity/ProppantPorosity.cpp
   )

set( dependencyList events dataRepository functions )

if( ENABLE_PVTPackage )
    set( constitutive_headers
         ${constitutive_headers}
         fluid/CompositionalMultiphaseFluid.hpp
         )

    set( constitutive_sources
         ${constitutive_sources}
         fluid/CompositionalMultiphaseFluid.cpp
         )
       
    add_subdirectory( PVTPackage )

    set( dependencyList ${dependencyList} PVTPackage )
endif()

if ( ENABLE_CUDA )
  set( dependencyList ${dependencyList} cuda )
endif()

blt_add_library( NAME       constitutive
                 SOURCES    ${constitutive_sources}
                 HEADERS    ${constitutive_headers}
                 DEPENDS_ON ${dependencyList}
                 OBJECT     ${GEOSX_BUILD_OBJ_LIBS}
               )
               
target_include_directories( constitutive PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )

add_subdirectory( unitTests )

geosx_add_code_checks( PREFIX constitutive 
                       EXCLUDES PVTPackage )
