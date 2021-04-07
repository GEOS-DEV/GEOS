#
# Specify all sources
#
set( constitutive_sources
     ConstitutiveManager.cpp
     ConstitutiveBase.cpp
     NullModel.cpp
     capillaryPressure/CapillaryPressureBase.cpp
     capillaryPressure/BrooksCoreyCapillaryPressure.cpp
     capillaryPressure/VanGenuchtenCapillaryPressure.cpp
     contact/ContactRelationBase.cpp
     contact/Coulomb.cpp
     fluid/DeadOilFluid.cpp
     fluid/CompressibleSinglePhaseFluid.cpp
     fluid/MultiPhaseMultiComponentFluid.cpp
     fluid/PVTFunctions/BrineCO2Density.cpp
     fluid/PVTFunctions/BrineViscosity.cpp
     fluid/PVTFunctions/CO2Solubility.cpp
     fluid/PVTFunctions/FenghourCO2Viscosity.cpp
     fluid/PVTFunctions/SpanWagnerCO2Density.cpp
     fluid/SingleFluidBase.cpp
     fluid/MultiFluidBase.cpp
     fluid/CompressibleSinglePhaseFluid.cpp
     fluid/SlurryFluidBase.cpp
     fluid/ProppantSlurryFluid.cpp
     fluid/ParticleFluidBase.cpp
     fluid/ParticleFluid.cpp               
     solid/Damage.cpp
     solid/DamageVolDev.cpp
     solid/DamageSpectral.cpp
     solid/DruckerPrager.cpp
     solid/DruckerPragerExtended.cpp
     solid/ElasticIsotropic.cpp
     solid/ElasticTransverseIsotropic.cpp
     solid/PoroElastic.cpp
     solid/PoreVolumeCompressibleSolid.cpp
     solid/SolidBase.cpp
     relativePermeability/RelativePermeabilityBase.cpp
     relativePermeability/BrooksCoreyRelativePermeability.cpp
     relativePermeability/BrooksCoreyBakerRelativePermeability.cpp
     relativePermeability/VanGenuchtenBakerRelativePermeability.cpp
     relativePermeability/TableRelativePermeability.cpp     
     capillaryPressure/CapillaryPressureBase.cpp
     capillaryPressure/BrooksCoreyCapillaryPressure.cpp
     capillaryPressure/VanGenuchtenCapillaryPressure.cpp
   )

if( BUILD_OBJ_LIBS)
  set( dependencyList dataRepository)
else()
  set( dependencyList common )
endif()



if( ENABLE_PVTPackage )
    set( constitutive_sources
         ${constitutive_sources}
         fluid/MultiFluidPVTPackageWrapper.cpp
         fluid/CompositionalMultiphaseFluid.cpp
         fluid/BlackOilFluid.cpp
       )
    set( dependencyList ${dependencyList} PVTPackage )
endif()

if ( ENABLE_OPENMP )
  set( dependencyList ${dependencyList} openmp )
endif()

if ( ENABLE_CUDA )
  set( dependencyList ${dependencyList} cuda )
endif()

blt_add_library( NAME                  constitutive
                 SOURCES               ${constitutive_sources}
                 DEPENDS_ON            ${dependencyList}
                 OBJECT                ${buildAsObj}
               )
               
target_include_directories( constitutive PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )

add_subdirectory( unitTests )

geosx_add_code_checks( PREFIX constitutive )
