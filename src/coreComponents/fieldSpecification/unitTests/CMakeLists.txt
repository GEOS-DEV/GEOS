# Specify list of tests
set( gtest_geosx_tests
     testFieldSpecificationsEnums.cpp
     testPermeabilitySpecification.cpp )

# Basic dependency list: fieldSpecification library plus gtest and parallel deps
set( dependencyList ${parallelDeps} fieldSpecification gtest )

# Add gtest C++ based tests
foreach(test ${gtest_geosx_tests})
    get_filename_component( test_name ${test} NAME_WE )
    blt_add_executable( NAME ${test_name}
                        SOURCES ${test}
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON ${dependencyList} )

    geos_add_test( NAME ${test_name}
                   COMMAND ${test_name} )

endforeach()
