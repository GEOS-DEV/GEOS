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

# Add pvtdriver xml based tests
blt_add_executable( NAME testPVT
                    SOURCES testPVT.cpp
                    OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                    DEPENDS_ON ${dependencyList} "-lz" ${tplDependencyList} )

# Add reactivefluiddriver xml based tests
blt_add_executable( NAME testReactiveFluid
                    SOURCES testReactiveFluid.cpp
                    OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                    DEPENDS_ON ${dependencyList} "-lz" ${tplDependencyList} )
