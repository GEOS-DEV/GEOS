#
# Specify all headers
#
set( constitutive_headers
     ConstitutiveManager.hpp
     ConstitutiveBase.hpp
     ExponentialRelation.hpp
     fluid/SingleFluidBase.hpp
     fluid/MultiFluidBase.hpp
     fluid/CompressibleSinglePhaseFluid.hpp
     fluid/SlurryFluidBase.hpp
     fluid/ProppantSlurryFluid.hpp
     fluid/ParticleFluidBase.hpp
     fluid/ParticleFluid.hpp          
     solid/LinearElasticAnisotropic.hpp
     solid/LinearElasticIsotropic.hpp
     solid/PoroElastic.hpp
     solid/PoreVolumeCompressibleSolid.hpp
     solid/SolidBase.hpp
     relativePermeability/RelativePermeabilityBase.hpp
     relativePermeability/BrooksCoreyRelativePermeability.hpp
     relativePermeability/BrooksCoreyBakerRelativePermeability.hpp
     relativePermeability/VanGenuchtenBakerRelativePermeability.hpp
     capillaryPressure/CapillaryPressureBase.hpp
     capillaryPressure/BrooksCoreyCapillaryPressure.hpp
     capillaryPressure/VanGenuchtenCapillaryPressure.hpp
     contact/ContactRelationBase.hpp
   )
#
# Specify all sources
#
set( constitutive_sources
     ConstitutiveManager.cpp
     ConstitutiveBase.cpp
     fluid/SingleFluidBase.cpp
     fluid/MultiFluidBase.cpp
     fluid/CompressibleSinglePhaseFluid.cpp
     fluid/SlurryFluidBase.cpp
     fluid/ProppantSlurryFluid.cpp
     fluid/ParticleFluidBase.cpp
     fluid/ParticleFluid.cpp               
     solid/LinearElasticAnisotropic.cpp
     solid/LinearElasticIsotropic.cpp
     solid/PoroElastic.cpp
     solid/PoreVolumeCompressibleSolid.cpp
     solid/SolidBase.cpp
     relativePermeability/RelativePermeabilityBase.cpp
     relativePermeability/BrooksCoreyRelativePermeability.cpp
     relativePermeability/BrooksCoreyBakerRelativePermeability.cpp
     relativePermeability/VanGenuchtenBakerRelativePermeability.cpp
     capillaryPressure/CapillaryPressureBase.cpp
     capillaryPressure/BrooksCoreyCapillaryPressure.cpp
     capillaryPressure/VanGenuchtenCapillaryPressure.cpp
     contact/ContactRelationBase.cpp
   )

if( BUILD_OBJ_LIBS)
  set( dependencyList dataRepository)
else()
  set( dependencyList common )
endif()



if( ENABLE_PVTPackage )
    set( constitutive_headers
         ${constitutive_headers}
         fluid/MultiFluidPVTPackageWrapper.hpp
         fluid/CompositionalMultiphaseFluid.hpp
         fluid/BlackOilFluid.hpp
       )
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
                 HEADERS               ${constitutive_headers}
                 DEPENDS_ON            ${dependencyList}
                 OBJECT                ${buildAsObj}
               )
               
target_include_directories( constitutive PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )

add_subdirectory( unitTests )

geosx_add_code_checks( PREFIX constitutive )
