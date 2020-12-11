set( fileIO_vtk_tests
     testVTK.cpp )

set( dependencyList dataRepository mesh managers vtk gmock gtest )

foreach(test ${fileIO_vtk_tests})
  get_filename_component( test_name ${test} NAME_WE )

  blt_add_executable( NAME ${test_name}
                      SOURCES ${test}
                      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON ${dependencyList} )

#  if ( ENABLE_MPI )
#    set( nranks 2 )
#    blt_add_test( NAME ${test_name}
#                  COMMAND ${test_name} -x ${nranks}
#                  NUM_MPI_TASKS ${nranks} )
#  else()
    blt_add_test( NAME ${test_name}
                  COMMAND ${test_name} )
#  endif()
endforeach()
