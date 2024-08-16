#
# Specify all headers
#
set( constitutiveDrivers_headers
     fluid/multifluid/PVTDriver.hpp
     fluid/multifluid/PVTDriverRunTest.hpp
     fluid/multifluid/reactive/ReactiveFluidDriver.hpp
     relativePermeability/RelpermDriver.hpp
     relativePermeability/RelpermDriverRunTest.hpp
     solid/TriaxialDriver.hpp
   )
#
# Specify all sources
#
set( constitutiveDrivers_sources
     fluid/multifluid/PVTDriver.cpp
     fluid/multifluid/blackOil/PVTDriverRunTestDeadOilFluid.cpp
     fluid/multifluid/CO2Brine/PVTDriverRunTestCO2BrinePhillipsFluid.cpp
     fluid/multifluid/CO2Brine/PVTDriverRunTestCO2BrinePhillipsThermalFluid.cpp
     fluid/multifluid/CO2Brine/PVTDriverRunTestCO2BrineEzrokhiFluid.cpp
     fluid/multifluid/CO2Brine/PVTDriverRunTestCO2BrineEzrokhiThermalFluid.cpp
     fluid/multifluid/compositional/PVTDriverRunTestCompositionalTwoPhaseConstantViscosity.cpp
     fluid/multifluid/compositional/PVTDriverRunTestCompositionalTwoPhaseLohrenzBrayClarkViscosity.cpp
     fluid/multifluid/reactive/ReactiveFluidDriver.cpp
     relativePermeability/RelpermDriver.cpp
     relativePermeability/RelpermDriverBrooksCoreyBakerRunTest.cpp
     relativePermeability/RelpermDriverBrooksCoreyStone2RunTest.cpp
     relativePermeability/RelpermDriverBrooksCoreyRunTest.cpp
     relativePermeability/RelpermDriverVanGenuchtenBakerRunTest.cpp
     relativePermeability/RelpermDriverVanGenuchtenStone2RunTest.cpp
     relativePermeability/RelpermDriverTableRelativeRunTest.cpp
     relativePermeability/RelpermDriverTableRelativeHysteresisRunTest.cpp
   )

set( dependencyList ${parallelDeps} constitutive events )

if( ENABLE_PVTPackage )
    set( constitutiveDrivers_sources
         ${constitutiveDrivers_sources}
          fluid/multifluid/compositional/PVTDriverRunTestCompositionalMultiphaseFluid.cpp
          )
endif()

blt_add_library( NAME       constitutiveDrivers
                 SOURCES    ${constitutiveDrivers_sources}
                 HEADERS    ${constitutiveDrivers_headers}
                 DEPENDS_ON ${dependencyList}
                 OBJECT     ${GEOS_BUILD_OBJ_LIBS}
                 SHARED     ${GEOS_BUILD_SHARED_LIBS}
               )

target_include_directories( constitutiveDrivers PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )

