#
# Specify list of tests
#

set(gtest_geosx_tests
    testDataTypes.cpp
    testTypeDispatch.cpp
    testFixedSizeDeque.cpp
   )

if ( ENABLE_CALIPER )
  list( APPEND gtest_geosx_tests
        testCaliperSmoke.cpp )
endif()

set( dependencyList common hdf5 gtest )

list( APPEND gtest_geosx_tests testLifoStorage.cpp )
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
