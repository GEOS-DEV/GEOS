#
# Specify all headers
#
set( constitutive_headers
     ConstitutiveManager.hpp
     ConstitutiveBase.hpp
     ExponentialRelation.hpp
     LinearElasticIsotropic.hpp
     Fluid/SingleFluidBase.hpp
     Fluid/MultiFluidBase.hpp
     Fluid/CompressibleSinglePhaseFluid.hpp
     Fluid/MultiPhaseMultiComponentFluid.hpp
     Fluid/PVTFunctions/BrineCO2DensityFunction.hpp
     Fluid/PVTFunctions/BrineViscosityFunction.hpp
     Fluid/PVTFunctions/CO2SolubilityFunction.hpp     
     Fluid/PVTFunctions/FenghourCO2ViscosityFunction.hpp
     Fluid/PVTFunctions/FlashModelBase.hpp
     Fluid/PVTFunctions/PVTFunctionBase.hpp
     Fluid/PVTFunctions/SpanWagnerCO2DensityFunction.hpp    
     Fluid/PVTFunctions/UtilityFunctions.hpp     
     Solid/PoreVolumeCompressibleSolid.hpp
     RelPerm/RelativePermeabilityBase.hpp
     RelPerm/BrooksCoreyRelativePermeability.hpp
     CapillaryPressure/CapillaryPressureBase.hpp
     CapillaryPressure/BrooksCoreyCapillaryPressure.hpp
     CapillaryPressure/VanGenuchtenCapillaryPressure.hpp
   )
#
# Specify all sources
#
set( constitutive_sources
     ConstitutiveManager.cpp
     ConstitutiveBase.cpp
     LinearElasticIsotropic.cpp
     Fluid/SingleFluidBase.cpp
     Fluid/MultiFluidBase.cpp
     Fluid/CompressibleSinglePhaseFluid.cpp
     Fluid/MultiPhaseMultiComponentFluid.cpp     
     Fluid/PVTFunctions/BrineCO2DensityFunction.cpp
     Fluid/PVTFunctions/BrineViscosityFunction.cpp 
     Fluid/PVTFunctions/CO2SolubilityFunction.cpp
     Fluid/PVTFunctions/FenghourCO2ViscosityFunction.cpp
     Fluid/PVTFunctions/FlashModelBase.cpp
     Fluid/PVTFunctions/PVTFunctionBase.cpp
     Fluid/PVTFunctions/SpanWagnerCO2DensityFunction.cpp    
     Fluid/PVTFunctions/UtilityFunctions.cpp
     Solid/PoreVolumeCompressibleSolid.cpp
     RelPerm/RelativePermeabilityBase.cpp
     RelPerm/BrooksCoreyRelativePermeability.cpp
     CapillaryPressure/CapillaryPressureBase.cpp
     CapillaryPressure/BrooksCoreyCapillaryPressure.cpp
     CapillaryPressure/VanGenuchtenCapillaryPressure.cpp
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


blt_add_library( NAME                  constitutive
                 SOURCES               ${constitutive_sources}
                 HEADERS               ${constitutive_headers}
                 DEPENDS_ON            ${dependencyList}
               )
               
target_include_directories( constitutive PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )

add_subdirectory( unitTests )

geosx_add_code_checks( PREFIX constitutive )
