#
# Specify all headers
#
set( constitutive_headers
     ConstitutiveManager.hpp
     ConstitutiveBase.hpp
     ExponentialRelation.hpp
     ExponentialRelation_impl.hpp
     LinearElasticIsotropic.hpp
     Fluid/MultiFluidBase.hpp
     Fluid/CompressibleSinglePhaseFluid.hpp
     Fluid/CompositionalMultiphaseFluid.hpp
     Solid/PoreVolumeCompressibleSolid.hpp
   )
#
# Specify all sources
#
set( constitutive_sources
     ConstitutiveManager.cpp
     ConstitutiveBase.cpp
     LinearElasticIsotropic.cpp
     Fluid/MultiFluidBase.cpp
     Fluid/CompressibleSinglePhaseFluid.cpp
     Fluid/CompositionalMultiphaseFluid.cpp
     Solid/PoreVolumeCompressibleSolid.cpp
   )

set( dependencyList dataRepository )

if(ENABLE_PVT_PACKAGE)
    set( dependencyList ${dependencyList} PVTPackage )
endif()


blt_add_library( NAME                  constitutive
                 SOURCES               ${constitutive_sources}
                 HEADERS               ${constitutive_headers}
                 DEPENDS_ON            ${dependencyList}
               )
               
target_include_directories( constitutive PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)

geosx_add_code_checks(PREFIX constitutive )