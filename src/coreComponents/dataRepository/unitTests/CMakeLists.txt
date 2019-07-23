#
# Specify list of tests
#

set( dataRepository_tests
     testDefaultValue.cpp
     testSidreBasic.cpp
     testSidreExtended.cpp
     testViewWrapper.cpp
   )

set( dependencyList dataRepository gtest cxx-utilities )

#
# Add gtest C++ based tests
#
foreach(test ${dataRepository_tests})
    get_filename_component( test_name ${test} NAME_WE )
    message("processing test_name ${test_name}")
    blt_add_executable( NAME ${test_name}
                        SOURCES ${test}
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON ${dependencyList}
                        )

    blt_add_test( NAME ${test_name} 
                  COMMAND ${test_name}
                  )
                  
endforeach()
