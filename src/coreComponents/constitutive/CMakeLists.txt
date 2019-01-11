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
     Solid/PoreVolumeCompressibleSolid.hpp
     RelPerm/RelativePermeabilityBase.hpp
     RelPerm/BrooksCoreyRelativePermeability.hpp
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
     Solid/PoreVolumeCompressibleSolid.cpp
     RelPerm/RelativePermeabilityBase.cpp
     RelPerm/BrooksCoreyRelativePermeability.cpp
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
