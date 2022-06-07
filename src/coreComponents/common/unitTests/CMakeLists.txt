#
# Specify list of tests
#

set(gtest_geosx_tests
    testDataTypes.cpp
    testTypeDispatch.cpp
   )

set( dependencyList ${deviceDeps} common hdf5 gtest )

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
