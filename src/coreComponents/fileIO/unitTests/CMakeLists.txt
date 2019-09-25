#
# Specify list of tests
#

set(testSources
    testLASIO.cpp
   )

set(dependencyList gtest fileIO geosx_core )

set( LAS_FILE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/las_example.las )
configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/wellFileNames.hpp.in ${CMAKE_BINARY_DIR}/include/tests/wellFileNames.hpp )

#
# Add gtest C++ based tests
#
foreach(test ${testSources})
    get_filename_component( test_name ${test} NAME_WE )
    blt_add_executable( NAME ${test_name}
                        SOURCES ${test}
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON ${dependencyList} )

    blt_add_test( NAME ${test_name} 
                  COMMAND ${test_name}
                  )

endforeach()
