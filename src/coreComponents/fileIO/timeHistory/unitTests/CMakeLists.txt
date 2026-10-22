# Specify list of tests
set( gtest_geosx_tests
     testHDFFile.cpp
     testHDFParallelFile.cpp )

set( dependencyList mainInterface gtest ${parallelDeps} HDF5::HDF5 )

# Add gtest C++ based tests
foreach(test ${gtest_geosx_tests})
    get_filename_component( test_name ${test} NAME_WE )
    blt_add_executable( NAME ${test_name}
                        SOURCES ${test}
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON ${dependencyList} )

    # Run parallel test with 2 ranks if MPI enabled
    if( ENABLE_MPI AND test_name STREQUAL "testHDFParallelFile" )
        geos_add_test( NAME ${test_name}
                       COMMAND ${test_name}
                       NUM_MPI_TASKS 2 )
    else()
        geos_add_test( NAME ${test_name}
                       COMMAND ${test_name} )
    endif()
endforeach()
