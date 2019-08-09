#
# Specify all headers
#
set( constitutive_headers
     ConstitutiveManager.hpp
     ConstitutiveBase.hpp
     ExponentialRelation.hpp
     Fluid/SingleFluidBase.hpp
     Fluid/MultiFluidBase.hpp
     Fluid/CompressibleSinglePhaseFluid.hpp
     Solid/LinearElasticAnisotropic.hpp
     Solid/LinearElasticIsotropic.hpp
     Solid/PoroElastic.hpp
     Solid/PoreVolumeCompressibleSolid.hpp
     Solid/SolidBase.hpp
     RelPerm/RelativePermeabilityBase.hpp
     RelPerm/BrooksCoreyRelativePermeability.hpp
     RelPerm/BrooksCoreyBakerRelativePermeability.hpp
     RelPerm/VanGenuchtenBakerRelativePermeability.hpp
     CapillaryPressure/CapillaryPressureBase.hpp
     CapillaryPressure/BrooksCoreyCapillaryPressure.hpp
     CapillaryPressure/VanGenuchtenCapillaryPressure.hpp
     contactRelations/ContactRelationBase.hpp
   )
#
# Specify all sources
#
set( constitutive_sources
     ConstitutiveManager.cpp
     ConstitutiveBase.cpp
     Fluid/SingleFluidBase.cpp
     Fluid/MultiFluidBase.cpp
     Fluid/CompressibleSinglePhaseFluid.cpp
     Solid/LinearElasticAnisotropic.cpp
     Solid/LinearElasticIsotropic.cpp
     Solid/PoroElastic.cpp
     Solid/PoreVolumeCompressibleSolid.cpp
     Solid/SolidBase.cpp
     RelPerm/RelativePermeabilityBase.cpp
     RelPerm/BrooksCoreyRelativePermeability.cpp
     RelPerm/BrooksCoreyBakerRelativePermeability.cpp
     RelPerm/VanGenuchtenBakerRelativePermeability.cpp
     CapillaryPressure/CapillaryPressureBase.cpp
     CapillaryPressure/BrooksCoreyCapillaryPressure.cpp
     CapillaryPressure/VanGenuchtenCapillaryPressure.cpp
     contactRelations/ContactRelationBase.cpp
   )

set( dependencyList dataRepository )

if( ENABLE_PVTPackage )
    set( constitutive_headers
         ${constitutive_headers}
         Fluid/MultiFluidPVTPackageWrapper.hpp
         Fluid/CompositionalMultiphaseFluid.hpp
         Fluid/BlackOilFluid.hpp
       )
    set( constitutive_sources
         ${constitutive_sources}
         Fluid/MultiFluidPVTPackageWrapper.cpp
         Fluid/CompositionalMultiphaseFluid.cpp
         Fluid/BlackOilFluid.cpp
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
                 OBJECT                TRUE
               )
               
target_include_directories( constitutive PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )

add_subdirectory( unitTests )

geosx_add_code_checks( PREFIX constitutive )
