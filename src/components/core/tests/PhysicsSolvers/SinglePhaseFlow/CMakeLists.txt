#
# Specify list of tests
#

set(gtest_geosx_tests
    testSinglePhaseFlow.cpp
    )

#
# Add gtest C++ based tests
#
foreach(test ${gtest_geosx_tests})
    get_filename_component( test_name ${test} NAME_WE )
    message("processing test_name ${test_name}")
    blt_add_executable( NAME ${test_name}
            SOURCES ${test}
            OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
            DEPENDS_ON ${geosx_core_libs}  gtest
            )

    blt_add_test( NAME ${test_name}
            COMMAND ${test_name}
            )

endforeach()