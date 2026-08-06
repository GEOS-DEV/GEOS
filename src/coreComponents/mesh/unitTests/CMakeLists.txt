# Specify list of tests
set( mesh_tests
     testMeshObjectPath.cpp
     testComputationalGeometry.cpp
     testGeometricObjects.cpp
     testGraphColoring.cpp
     testMeshEnums.cpp )

if( ENABLE_METIS )
  list( APPEND mesh_tests testMETISInterface.cpp )
endif()

set( mesh_tests_mpi
     testGraphColoringMPI.cpp )



set( dependencyList blas lapack gtest mesh ${parallelDeps} )

# Add gtest C++ based tests
foreach(test ${mesh_tests})
    get_filename_component( test_name ${test} NAME_WE )
    blt_add_executable( NAME ${test_name}
                        SOURCES ${test}
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON ${dependencyList} )

    geos_add_test( NAME ${test_name}
                   COMMAND ${test_name} )
endforeach()






if( ENABLE_MPI )
set( dependencyList blas lapack gtest mesh ${parallelDeps} )

if (ENABLE_TRILINOS)
   list( APPEND dependencyList trilinos )
endif()

set( nranks 4 )

  foreach( test ${mesh_tests_mpi} )
    get_filename_component( file_we ${test} NAME_WE )
    set( test_name ${file_we}_mpi )
    blt_add_executable( NAME ${test_name}
                        SOURCES ${test}
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON ${dependencyList} ${tplDependencyList} )

    geos_add_test( NAME ${test_name}
                   COMMAND ${test_name} -x ${nranks}
                   NUM_MPI_TASKS ${nranks} )
    # Older g++ version need a hint to use std::filesystem
    if( CMAKE_CXX_COMPILER_ID STREQUAL "GNU" )
        if( CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9 )
            target_link_libraries( ${test_name} PRIVATE "stdc++fs" )
        endif ()
    endif ()
    if ( ENABLE_CUDA )
      target_link_libraries( ${test_name} PRIVATE "stdc++fs" )
    endif ()
  endforeach()

  # VTKMeshScattering test — 4 MPI ranks
  blt_add_executable( NAME testVTKMeshScattering_mpi
                      SOURCES testVTKMeshScattering.cpp
                      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON ${dependencyList} ${tplDependencyList} )

  geos_add_test( NAME testVTKMeshScattering_mpi
                 COMMAND testVTKMeshScattering_mpi
                 NUM_MPI_TASKS 4 )
endif()
