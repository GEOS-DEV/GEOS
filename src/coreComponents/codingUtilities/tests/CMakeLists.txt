#
# Specify list of tests
#

set( testSources
     testGeosxTraits.cpp
     testParallelTestUtilities.cpp
     )

set( testParallelSources
     testParallelTestUtilities.cpp
     )

set( dependencyList gtest codingUtilities )

if( GEOSX_BUILD_SHARED_LIBS )
  set( dependencyList ${dependencyList} geosx_core )
else()
  set( dependencyList ${dependencyList} ${geosx_core_libs} )
endif()

if( ENABLE_CUDA )
  set( dependencyList ${dependencyList} cuda )
endif()

#
# Add gtest C++ based tests
#
foreach( test ${testSources} )
  get_filename_component( test_name ${test} NAME_WE )
  blt_add_executable( NAME ${test_name}
                      SOURCES ${test}
                      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON ${dependencyList} )

  blt_add_test( NAME ${test_name}
                COMMAND ${test_name}
                )
endforeach()

if( ENABLE_MPI )

  set( nranks 2 )

  foreach( test ${testParallelSources} )
    get_filename_component( file_we ${test} NAME_WE )
    set( test_name ${file_we}_mpi )
    blt_add_executable( NAME ${test_name}
                        SOURCES ${test}
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON ${dependencyList}
                        )

    blt_add_test( NAME ${test_name}
                  COMMAND ${test_name} -x ${nranks}
                  NUM_MPI_TASKS ${nranks}
                  )
  endforeach()
endif()
