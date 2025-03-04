# Specify list of tests
set( gtest_geosx_tests
     testCapillaryPressure.cpp
     testCO2BrinePVTModels.cpp
     testCO2SpycherPruessModels.cpp
     testDamage.cpp
     testMultiFluidCO2Brine.cpp
     testMultiFluidCompositionalMultiphase.cpp
     testMultiFluidDeadOil.cpp
     testMultiFluidLiveOil.cpp
     testRelPerm.cpp
     testRelPermHysteresis.cpp )

set( gtest_triaxial_xmls
     driverTests/testTriaxial_druckerPragerExtended.xml
     driverTests/testTriaxial_elasticIsotropic.xml
     driverTests/testTriaxial_elasticIsotropicPressureDependent.xml
     driverTests/testTriaxial_modifiedCamClay.xml
     driverTests/testTriaxial_delftEggUseLinear.xml
     driverTests/testTriaxial_modifiedCamClayVolumetric.xml
     driverTests/testTriaxial_delftEggLoadPathDryUseLinear.xml
     driverTests/testTriaxial_delftEggLoadPathWetUseLinear.xml
     driverTests/testTriaxial_delftEggCase1.xml
     driverTests/testTriaxial_delftEggCase2.xml )

set( gtest_pvt_xmls
     driverTests/testPVT_docExample.xml
     driverTests/testPVT.xml
     driverTests/testPVT_CO2Brine.xml
     driverTests/testPVT_CO2BrineTables.xml
     driverTests/testPVT_PhaseComposition.xml
     driverTests/testPVT_ThreePhaseCompositional.xml
   )

set( gtest_reactivefluid_xmls
     driverTests/testReactiveFluid.xml )

set( tplDependencyList ${parallelDeps} gtest )

set( dependencyList mainInterface  )

if( ENABLE_PVTPackage )
    list( APPEND gtest_geosx_tests
          testMultiFluidCompositionalMultiphasePVTPackage.cpp )
endif()

geos_decorate_link_dependencies( LIST decoratedDependencies
                                 DEPENDENCIES ${dependencyList} )

# Add gtest C++ based tests
foreach(test ${gtest_geosx_tests})
    get_filename_component( test_name ${test} NAME_WE )
    blt_add_executable( NAME ${test_name}
                        SOURCES ${test}
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON ${decoratedDependencies} ${tplDependencyList} )
    geos_add_test( NAME ${test_name}
                   COMMAND ${test_name} )

endforeach()

# Add triaxial xml based tests
blt_add_executable( NAME testTriaxial
                    SOURCES testTriaxial.cpp
                    OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                    DEPENDS_ON ${dependencyList} "-lz" ${tplDependencyList})

foreach(test ${gtest_triaxial_xmls})
    get_filename_component( test_name ${test} NAME_WE )
    geos_add_test( NAME ${test_name}
                   COMMAND testTriaxial -i ${CMAKE_CURRENT_LIST_DIR}/${test} )
endforeach()

# Add pvtdriver xml based tests
blt_add_executable( NAME testPVT
                    SOURCES testPVT.cpp
                    OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                    DEPENDS_ON ${dependencyList} "-lz" ${tplDependencyList} )

foreach(test ${gtest_pvt_xmls})
    get_filename_component( test_name ${test} NAME_WE )
    geos_add_test( NAME ${test_name}
                   COMMAND testPVT -i ${CMAKE_CURRENT_LIST_DIR}/${test} )
endforeach()


# Add reactivefluiddriver xml based tests
blt_add_executable( NAME testReactiveFluid
                    SOURCES testReactiveFluid.cpp
                    OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                    DEPENDS_ON ${dependencyList} "-lz" ${tplDependencyList} )

foreach(test ${gtest_reactivefluid_xmls})
    get_filename_component( test_name ${test} NAME_WE )
    geos_add_test( NAME ${test_name}
                   COMMAND testReactiveFluid -i ${CMAKE_CURRENT_LIST_DIR}/${test} )
endforeach()