#
# Specify list of tests
#


set( gtest_geosx_tests
     testMeshEnums.cpp
     testMeshGeneration.cpp
     testNeighborCommunicator.cpp
     )

set( gtest_geosx_mpi_tests
     testNeighborCommunicator.cpp
     )

if( ENABLE_PAMELA )
  list( APPEND gtest_geosx_tests
        testPAMELAImport.cpp
        )
endif()

if( ENABLE_VTK )
  list( APPEND gtest_geosx_tests
        testVTKImport.cpp
        )
  list( APPEND gtest_geosx_mpi_tests
        testVTKImport.cpp
        )
endif()

if( ENABLE_VTK OR ENABLE_PAMELA )
  set( TEST_MESH_DIR ${CMAKE_CURRENT_SOURCE_DIR} )
  configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/meshDirName.hpp.in ${CMAKE_BINARY_DIR}/include/tests/meshDirName.hpp )
endif()

set( dependencyList gtest )

if( GEOSX_BUILD_SHARED_LIBS )
  set( dependencyList ${dependencyList} geosx_core )
else()
  set( dependencyList ${dependencyList} ${geosx_core_libs} )
endif()

if( ENABLE_CUDA )
  set( dependencyList ${dependencyList} cuda )
endif()

if( ENABLE_CUDA_NVTOOLSEXT )
  set( dependencyList ${dependencyList} CUDA::nvToolsExt )
endif()

if( ENABLE_OPENMP )
  set( dependencyList ${dependencyList} openmp )
endif()

#
# Add gtest C++ based tests
#
foreach( test ${gtest_geosx_tests} )
  get_filename_component( test_name ${test} NAME_WE )
  blt_add_executable( NAME ${test_name}
                      SOURCES ${test}
                      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON ${dependencyList}
                      )

  blt_add_test( NAME ${test_name}
                COMMAND ${test_name} ${CMAKE_CURRENT_LIST_DIR}
                )
endforeach()

if( ENABLE_MPI )

  set( nranks 2 )

  foreach( test ${gtest_geosx_mpi_tests} )
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
