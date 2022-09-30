#
# Specify list of tests
#

set( gtest_geosx_tests
     testDamage.cpp
     testRelPerm.cpp
     testRelPermHysteresis.cpp
     testCapillaryPressure.cpp
   )

set( gtest_triaxial_xmls
     testTriaxial_druckerPragerExtended.xml
     testTriaxial_elasticIsotropic.xml
     testTriaxial_elasticIsotropicPressureDependent.xml
     testTriaxial_modifiedCamClay.xml
     testTriaxial_delftEggUseLinear.xml
     testTriaxial_modifiedCamClayVolumetric.xml
     testTriaxial_delftEggLoadPathDryUseLinear.xml
     testTriaxial_delftEggLoadPathWetUseLinear.xml
     testTriaxial_delftEggCase1.xml
     testTriaxial_delftEggCase2.xml
   )

set( gtest_pvt_xmls
     testPVT.xml
   )

set( input_data
     testPVT_brine.txt
     testPVT_hydrogenMixtureB.txt
     testTriaxial_delftEggCase2.txt
     testTriaxial_delftEggLoadPathWetUseLinear.txt
     testTriaxial_druckerPragerExtended.txt
     testTriaxial_elasticIsotropicPressureDependent.txt
     testTriaxial_modifiedCamClayVolumetric.txt
     testPVT_CO2.txt
     testPVT_hydrogenMixtureA.txt
     testTriaxial_delftEggCase1.txt
     testTriaxial_delftEggLoadPathDryUseLinear.txt
     testTriaxial_delftEggUseLinear.txt
     testTriaxial_elasticIsotropic.txt
     testTriaxial_modifiedCamClay.txt
     testPVT_data/brinePVT.txt
     testPVT_data/oilPVT.txt
     testPVT_data/carbonDioxideFlash.txt
     testPVT_data/referencePVT_CO2.txt
     testPVT_data/carbonDioxidePVT.txt
     testPVT_data/referencePVT_H2_CH4.txt
     testPVT_data/gasPVT.txt
     testPVT_data/referencePVT_brine.txt
  )

foreach( file ${gtest_triaxial_xmls} ${gtest_pvt_xmls} ${input_data} )
  configure_file( ${file} ${CMAKE_CURRENT_BINARY_DIR}/${file} COPYONLY )
endforeach()

set( dependencyList gtest )

if ( GEOSX_BUILD_SHARED_LIBS )
  set (dependencyList ${dependencyList} geosx_core )
else()
  set (dependencyList ${dependencyList} ${geosx_core_libs} )
endif()

if( ENABLE_PVTPackage )
    list( APPEND gtest_geosx_tests
          testMultiFluid.cpp
          testCO2BrinePVTModels.cpp )

    set( dependencyList ${dependencyList} PVTPackage )
endif()

if( ENABLE_CUDA )
  set( dependencyList ${dependencyList} cuda )
endif()

if( ENABLE_CUDA_NVTOOLSEXT )
  set( dependencyList ${dependencyList} CUDA::nvToolsExt )
endif()

if( ENABLE_OPENMP )
  set( dependencyList ${dependencyList} openmp )
endif()

#
# Add gtest C++ based tests
#
foreach(test ${gtest_geosx_tests})
    get_filename_component( test_name ${test} NAME_WE )
    blt_add_executable( NAME ${test_name}
            SOURCES ${test}
            OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
            DEPENDS_ON ${dependencyList}
            )

    blt_add_test( NAME ${test_name}
                  COMMAND ${test_name}
                )

endforeach()

#
# Add triaxial xml based tests
#

blt_add_executable( NAME testTriaxial
                    SOURCES testTriaxial.cpp
                    OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                    DEPENDS_ON ${dependencyList} "-lz"
                  )

foreach(test ${gtest_triaxial_xmls})
    get_filename_component( test_name ${test} NAME_WE )
    blt_add_test( NAME ${test_name}
                  COMMAND testTriaxial -i ${CMAKE_CURRENT_BINARY_DIR}/${test}
                  )
endforeach()

#
# Add pvtdriver xml based tests
#

blt_add_executable( NAME testPVT
                    SOURCES testPVT.cpp
                    OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                    DEPENDS_ON ${dependencyList} "-lz"
                  )

foreach(test ${gtest_pvt_xmls})
    get_filename_component( test_name ${test} NAME_WE )
    blt_add_test( NAME ${test_name}
                  COMMAND testPVT -i ${CMAKE_CURRENT_BINARY_DIR}/${test}
                  )
endforeach()
