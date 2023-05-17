#
# Specify list of tests
#

set( gtest_geosx_tests
     testAquiferBoundaryCondition.cpp
     testFieldSpecificationsEnums.cpp
     testRecursiveFieldApplication.cpp
   )


set( dependencyList ${parallelDeps} gtest )

if ( GEOSX_BUILD_SHARED_LIBS )
  set (dependencyList ${dependencyList} geosx_core )
else()
  set (dependencyList ${dependencyList} ${geosx_core_libs} )
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
