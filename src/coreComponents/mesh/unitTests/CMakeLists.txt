#
# Specify list of tests
#

set( gtest_geosx_tests
     testGraphGeneration.cpp
   )

set( dependencyList gtest )

if ( GEOSX_BUILD_SHARED_LIBS )
  set (dependencyList ${dependencyList} geosx_core)
else()
  set (dependencyList ${dependencyList} ${geosx_core_libs} )
endif()

if ( ENABLE_MPI )
  set ( dependencyList ${dependencyList} mpi )
endif()

if( ENABLE_OPENMP )
  set( dependencyList ${dependencyList} openmp )
endif()

if ( ENABLE_CUDA )
  set( dependencyList ${dependencyList} cuda )
endif()
set( TXT_FILE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/BaseGraph.txt )
configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/graphFileNames.hpp.in ${CMAKE_BINARY_DIR}/include/tests/graphFileNames.hpp )


#
# Add gtest C++ based tests
#
foreach(test ${gtest_geosx_tests})
  get_filename_component( test_name ${test} NAME_WE )

  blt_add_executable( NAME ${test_name}
                      SOURCES ${test}
                      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON ${dependencyList} )

  blt_add_test( NAME ${test_name}
                COMMAND ${test_name} )
endforeach()

# For some reason, BLT is not setting CUDA language for these source files
if ( ENABLE_CUDA )
  set_source_files_properties( ${gtest_geosx_tests} PROPERTIES LANGUAGE CUDA )
endif()
