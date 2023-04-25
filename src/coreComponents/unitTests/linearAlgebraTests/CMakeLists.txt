# Specify list of tests
set( LAI_tests
     testDofManager.cpp
     testLAIHelperFunctions.cpp )

set( nranks 2 )

# Add gtest C++ based tests
set( dependencyList ${parallelDeps} gtest )

if ( GEOSX_BUILD_SHARED_LIBS )
  list( APPEND dependencyList geosx_core )
else()
  list( APPEND dependencyList ${geosx_core_libs} )
endif()

foreach(test ${LAI_tests})
  get_filename_component( test_name ${test} NAME_WE )

  blt_add_executable( NAME ${test_name}
                      SOURCES ${test}
                      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON ${dependencyList} pugixml::pugixml pugixml fmt::fmt
                    )

if ( ENABLE_MPI )
  geos_add_test( NAME ${test_name}
                 COMMAND ${test_name} -x ${nranks}
                 NUM_MPI_TASKS ${nranks} )
else()
  geos_add_test( NAME ${test_name}
                 COMMAND ${test_name} )
endif()
endforeach()
